// MIT License
//
// Copyright (c) 2025 Jussi Lind
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <algorithm>
#include <cmath>
#include <fftw3.h>
#include <iostream>
#include <numeric>
#include <sndfile.h>
#include <stdexcept>
#include <string>
#include <vector>

// --- Helper Structs and Functions ---

struct Band
{
    double low;
    double high;
    double f_center; // Geometric mean of low and high
};

// Map frequency in Hz to FFT bin index
size_t freqToBin(double freq, double fs, size_t N)
{
    if (N == 0 || fs == 0)
        return 0;
    return static_cast<size_t>(std::round(freq / fs * N));
}

// Calculates the geometric center frequency
double getBandCenterFrequency(double low, double high)
{
    return std::sqrt(low * high);
}

// Generates logarithmic bands based on the requested count
std::vector<Band> generate_bands(double samplerate, int num_bands)
{
    std::vector<Band> bands;

    if (num_bands < 1) {
        std::cerr << "Error: Number of bands must be at least 1.\n";
        return bands;
    }

    const double f_start = 20.0;
    double nyquist = samplerate / 2.0;
    double f_end = std::min(20000.0, nyquist - 1.0);

    if (f_end <= f_start) {
        std::cerr << "Warning: Samplerate (" << samplerate
                  << " Hz) is too low to define full spectral range (Nyquist=" << nyquist
                  << " Hz).\n";
        if (nyquist > 20.0) {
            bands.push_back({f_start, nyquist, getBandCenterFrequency(f_start, nyquist)});
        }
        return bands;
    }

    double multiplier = std::pow(f_end / f_start, 1.0 / num_bands);
    double low = f_start;

    for (int i = 0; i < num_bands; ++i) {
        double high = low * multiplier;
        if (i == num_bands - 1) {
            high = f_end; // Ensure the last band hits the ceiling
        }
        double f_center = getBandCenterFrequency(low, high);
        bands.push_back({low, high, f_center});
        low = high;
    }
    return bands;
}

// Reads a file, calculates its mono average, and stores full buffer data.
std::vector<double> read_file_and_get_mono_data(const std::string &filepath,
                                                SF_INFO &sfinfo_out,
                                                std::vector<double> &full_buffer_out)
{
    SNDFILE *inFile = sf_open(filepath.c_str(), SFM_READ, &sfinfo_out);
    if (!inFile) {
        throw std::runtime_error("Error opening file " + filepath + ": " + sf_strerror(NULL));
    }

    if (sfinfo_out.frames <= 0 || sfinfo_out.channels <= 0 || sfinfo_out.samplerate <= 0) {
        sf_close(inFile);
        throw std::runtime_error("Error: Invalid file info (frames, channels, or samplerate) for "
                                 + filepath + ".");
    }

    // Read all samples
    full_buffer_out.resize(sfinfo_out.frames * sfinfo_out.channels);
    sf_readf_double(inFile, full_buffer_out.data(), sfinfo_out.frames);
    sf_close(inFile);

    const sf_count_t N = sfinfo_out.frames;
    const int num_channels = sfinfo_out.channels;

    // Generate Mono Average Signal
    std::vector<double> mono_data(N);
    for (sf_count_t i = 0; i < N; ++i) {
        double sum = 0;
        for (int c = 0; c < num_channels; ++c) {
            sum += full_buffer_out[i * num_channels + c];
        }
        mono_data[i] = sum / num_channels;
    }
    return mono_data;
}

// Performs FFT on mono data and computes N-band RMS energy.
std::vector<double> calculate_band_amp(const std::vector<double> &mono_data,
                                       const SF_INFO &sfinfo_in,
                                       const std::vector<Band> &bands)
{
    const sf_count_t N = mono_data.size();

    if (N < 2) {
        throw std::runtime_error("Error: Sample size is too small for FFT.");
    }

    fftw_complex *fft_mono = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1));
    fftw_plan forward_mono = fftw_plan_dft_r2c_1d(N,
                                                  const_cast<double *>(mono_data.data()),
                                                  fft_mono,
                                                  FFTW_ESTIMATE);

    fftw_execute(forward_mono);

    // Compute Band Amplitude (proportional to RMS)
    std::vector<double> band_amp;
    const double fs = sfinfo_in.samplerate;

    for (const auto &band : bands) {
        size_t start = freqToBin(band.low, fs, N);
        size_t end = freqToBin(band.high, fs, N);

        if (end >= N / 2)
            end = N / 2;
        if (start >= end) {
            band_amp.push_back(1e-12); // Tiny floor for truly empty bands
            continue;
        }

        // Exclude DC bin (k=0) unless the band starts at 0 Hz
        if (start == 0 && band.low > 0)
            start = 1;

        double sum_sq_mag = 0;
        for (size_t k = start; k <= end; ++k) {
            double mag_sq = fft_mono[k][0] * fft_mono[k][0] + fft_mono[k][1] * fft_mono[k][1];
            // Correct scaling: bins 1 to N/2 - 1 are two-sided, so their energy is doubled.
            sum_sq_mag += (k > 0 && k < N / 2) ? 2.0 * mag_sq : mag_sq;
        }

        band_amp.push_back(
            std::sqrt(sum_sq_mag)); // No floor here, let the dynamic floor handle it later
    }

    fftw_destroy_plan(forward_mono);
    fftw_free(fft_mono);

    return band_amp;
}

int main(int argc, char **argv)
{
    // Usage: fft_match input.wav reference.wav output.wav [num_bands (def 8)] [max_boost_db (def 3.0)] [max_cut_db (def 3.0)]
    if (argc < 4 || argc > 7) {
        std::cout << "Usage: fft_match input.wav reference.wav output.wav [num_bands (default 8)] "
                     "[max_boost_db (default 3.0)] [max_cut_db (default 3.0)]\n";
        std::cout
            << "This tool balances the N-band spectral RMS of input.wav to match reference.wav.\n";
        return 1;
    }

    const std::string infile = argv[1];
    const std::string reffile = argv[2];
    const std::string outfile = argv[3];

    // --- Argument Parsing for Optional Parameters ---
    int num_bands_target = 8;
    double max_boost_db = 3.0; // Default max boost
    double max_cut_db = 3.0;   // Default max cut

    try {
        if (argc >= 5) {
            num_bands_target = std::stoi(argv[4]);
            if (num_bands_target < 2 || num_bands_target > 30) {
                std::cerr
                    << "Warning: Number of bands should be between 2 and 30. Using default 8.\n";
                num_bands_target = 8;
            }
        }

        if (argc >= 6) {
            max_boost_db = std::stod(argv[5]);
            if (max_boost_db < 0.0) {
                std::cerr << "Warning: Max boost dB cannot be negative. Using default 3.0.\n";
                max_boost_db = 3.0;
            }
        }

        if (argc == 7) {
            max_cut_db = std::stod(argv[6]);
            if (max_cut_db < 0.0) {
                std::cerr << "Warning: Max cut dB cannot be negative. Using default 3.0.\n";
                max_cut_db = 3.0;
            }
        }

    } catch (const std::exception &e) {
        std::cerr << "Warning: Invalid numeric argument provided for optional parameters. Using "
                     "defaults.\n";
        num_bands_target = 8;
        max_boost_db = 3.0;
        max_cut_db = 3.0;
    }
    // --- End Argument Parsing ---

    SF_INFO input_sfinfo = {};
    std::vector<double> input_buffer;
    std::vector<double> input_mono_data;

    SF_INFO ref_sfinfo = {};
    std::vector<double> ref_buffer;
    std::vector<double> ref_mono_data;

    std::vector<Band> bands;
    std::vector<double> input_band_amp;
    std::vector<double> ref_band_amp;

    try {
        // --- 1. Load Data and File Info (Input and Reference) ---
        std::cout << "Loading input track: " << infile << "...\n";
        input_mono_data = read_file_and_get_mono_data(infile, input_sfinfo, input_buffer);

        std::cout << "Loading reference track: " << reffile << "...\n";
        ref_mono_data = read_file_and_get_mono_data(reffile, ref_sfinfo, ref_buffer);

        if (input_sfinfo.samplerate != ref_sfinfo.samplerate) {
            std::cerr << "Error: Input and Reference files must have the same sample rate!\n";
            return 1;
        }

        // --- 2. Generate Bands ---
        bands = generate_bands(input_sfinfo.samplerate, num_bands_target);

        if (bands.empty()) {
            std::cerr
                << "Error: Failed to generate valid frequency bands for the given sample rate ("
                << input_sfinfo.samplerate << " Hz).\n";
            return 1;
        }

        const int num_channels = input_sfinfo.channels;
        const sf_count_t num_frames = input_sfinfo.frames;
        const size_t N = num_frames;

        std::cout << "Spectral analysis parameters: Sample Rate=" << input_sfinfo.samplerate
                  << " Hz, Channels=" << num_channels << ", Frames=" << num_frames
                  << ", Bands=" << bands.size() << "\n";

        // --- 3. Calculate Amplitudes for both tracks ---
        input_band_amp = calculate_band_amp(input_mono_data, input_sfinfo, bands);
        ref_band_amp = calculate_band_amp(ref_mono_data, ref_sfinfo, bands);

        std::cout << "Spectral analysis complete. Calculating gains...\n";

        // --- 4. GAIN CALCULATION (WITH LEVEL NORMALIZATION & DYNAMIC FLOOR) ---

        // Calculate overall spectral loudness (average of band amplitudes)
        double input_overall_amp = std::accumulate(input_band_amp.begin(), input_band_amp.end(), 0.0)
                                   / input_band_amp.size();
        double ref_overall_amp = std::accumulate(ref_band_amp.begin(), ref_band_amp.end(), 0.0)
                                 / ref_band_amp.size();

        // Dynamic Amplitude Floor (e.g., -70 dB relative to the louder overall average band RMS)
        // 1e-7 is -70dB (20 * log10(1e-7) = -140dB, but since we compare RMS values, which are magnitude not power,
        // 1e-7 is 7 orders of magnitude down in amplitude.)
        const double overall_max_amp = std::max(input_overall_amp, ref_overall_amp);
        const double AMPLITUDE_FLOOR = std::max(overall_max_amp * 1e-7, 1e-12);

        // Level match factor
        const double level_match_factor = (ref_overall_amp > 1e-12)
                                              ? input_overall_amp / ref_overall_amp
                                              : 1.0;

        std::cout << "Input Overall Avg Band RMS: " << input_overall_amp << "\n";
        std::cout << "Reference Overall Avg Band RMS: " << ref_overall_amp << "\n";
        std::cout << "Applying level match factor (Ref * " << level_match_factor
                  << ") to normalize spectra before EQ.\n";
        std::cout << "Using Dynamic Amplitude Floor (to prevent silence errors): "
                  << AMPLITUDE_FLOOR << "\n";

        std::vector<double> gains;

        // Calculate the absolute gain factors based on user-defined dB limits
        const double MAX_GAIN = std::pow(10.0, max_boost_db / 20.0);
        const double MIN_GAIN = std::pow(10.0, -max_cut_db / 20.0);

        std::cout << "Using Gain Limits: Boost=" << max_boost_db << " dB (factor " << MAX_GAIN
                  << "), Cut=" << max_cut_db << " dB (factor " << MIN_GAIN << ").\n";

        for (size_t b = 0; b < bands.size(); ++b) {
            double target_amp_raw = ref_band_amp[b] * level_match_factor;

            // Apply the dynamic floor for stability before calculating gain
            double current_amp_floored = std::max(input_band_amp[b], AMPLITUDE_FLOOR);
            double target_amp_floored = std::max(target_amp_raw, AMPLITUDE_FLOOR);

            // Gain is Target / Current
            double G = target_amp_floored / current_amp_floored;

            // Apply hard limits
            if (G > MAX_GAIN)
                G = MAX_GAIN; // Limits boost
            if (G < MIN_GAIN)
                G = MIN_GAIN; // Limits cut

            gains.push_back(G);

            double dB = 20.0 * std::log10(G);
            std::cout << "Band " << b + 1 << " (" << bands[b].low << " Hz - " << bands[b].high
                      << " Hz): Gain = " << dB << " dB\n";
        }

        // --- 5. PROCESS EACH CHANNEL (FFT and IFFT application remains the same) ---

        std::vector<std::vector<double>> processed_channels(num_channels,
                                                            std::vector<double>(num_frames));
        double max_val_global = 0.0;

        std::vector<double> channel_data(num_frames);
        fftw_complex *fft_channel = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1));
        fftw_complex *ifft_channel = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * N);
        fftw_plan forward_channel = fftw_plan_dft_r2c_1d(N,
                                                         channel_data.data(),
                                                         fft_channel,
                                                         FFTW_ESTIMATE);
        fftw_plan backward_channel = fftw_plan_dft_c2r_1d(N,
                                                          ifft_channel,
                                                          channel_data.data(),
                                                          FFTW_ESTIMATE);

        std::cout << "Applying EQ to " << num_channels << " channels...\n";

        for (int c = 0; c < num_channels; ++c) {
            // Extract channel data
            for (sf_count_t i = 0; i < num_frames; ++i) {
                channel_data[i] = input_buffer[i * num_channels + c];
            }

            // Forward FFT
            fftw_execute(forward_channel);

            // Copy spectrum to IFFT buffer
            for (size_t i = 0; i <= N / 2; i++) {
                ifft_channel[i][0] = fft_channel[i][0];
                ifft_channel[i][1] = fft_channel[i][1];
            }

            // Apply gains to the positive frequency side
            for (size_t b = 0; b < bands.size(); ++b) {
                double G = gains[b];
                size_t start = freqToBin(bands[b].low, input_sfinfo.samplerate, N);
                size_t end = freqToBin(bands[b].high, input_sfinfo.samplerate, N);

                if (end >= N / 2)
                    end = N / 2;
                if (start >= end)
                    continue;

                if (start == 0 && bands[b].low > 0)
                    start = 1;

                for (size_t k = start; k <= end; ++k) {
                    ifft_channel[k][0] *= G;
                    ifft_channel[k][1] *= G;
                }
            }

            // Set the negative frequency side (symmetry)
            for (size_t k = 1; k < N / 2; ++k) {
                size_t nk = N - k;
                ifft_channel[nk][0] = ifft_channel[k][0];
                ifft_channel[nk][1] = -ifft_channel[k][1];
            }

            // Inverse FFT
            fftw_execute(backward_channel);

            // Normalize (FFTW inverse scales by N) and track global max
            for (sf_count_t i = 0; i < num_frames; ++i) {
                double v = channel_data[i] / N;
                processed_channels[c][i] = v;
                if (std::abs(v) > max_val_global)
                    max_val_global = std::abs(v);
            }
        }

        // Cleanup Channel FFTW resources
        fftw_destroy_plan(forward_channel);
        fftw_destroy_plan(backward_channel);
        fftw_free(fft_channel);
        fftw_free(ifft_channel);

        // --- 6. Final Interleave and Write ---

        if (max_val_global < 1e-12)
            max_val_global = 1.0;

        // Apply final normalization
        const double final_scale = (max_val_global > 1.0) ? 1.0 / max_val_global : 1.0;
        std::cout << "Applying final output scaling factor: " << final_scale << "\n";

        std::vector<double> final_interleaved_data(num_frames * num_channels);

        for (sf_count_t i = 0; i < num_frames; ++i) {
            for (int c = 0; c < num_channels; ++c) {
                final_interleaved_data[i * num_channels + c] = processed_channels[c][i]
                                                               * final_scale;
            }
        }

        // Write multi-channel output
        SF_INFO outInfo = input_sfinfo;
        outInfo.channels = num_channels;

        SNDFILE *outFile = sf_open(outfile.c_str(), SFM_WRITE, &outInfo);
        if (!outFile) {
            std::cerr << "Error opening output file: " << sf_strerror(NULL) << "\n";
            return 1;
        }

        sf_writef_double(outFile, final_interleaved_data.data(), num_frames);
        sf_close(outFile);

        std::cout << "Successfully matched spectral balance and wrote to: " << outfile << "\n";
        return 0;

    } catch (const std::runtime_error &e) {
        std::cerr << e.what() << "\n";
        return 1;
    } catch (...) {
        std::cerr << "An unknown error occurred.\n";
        return 1;
    }
}
