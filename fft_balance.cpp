#include "fft_balance.hpp"
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <numeric>
#include <algorithm>

FFTBalance::FFTBalance(int num_bands, double max_boost_db, double max_cut_db)
    : m_num_bands_target(num_bands), m_max_boost_db(max_boost_db), m_max_cut_db(max_cut_db),
      m_forward_channel_plan(nullptr), m_backward_channel_plan(nullptr),
      m_fft_channel(nullptr), m_ifft_channel(nullptr)
{
}

FFTBalance::~FFTBalance()
{
    if (m_forward_channel_plan) fftw_destroy_plan(m_forward_channel_plan);
    if (m_backward_channel_plan) fftw_destroy_plan(m_backward_channel_plan);
    if (m_fft_channel) fftw_free(m_fft_channel);
    if (m_ifft_channel) fftw_free(m_ifft_channel);
}

void FFTBalance::load_tracks(const std::string& infile, const std::string& reffile)
{
    std::cout << "Loading input track: " << infile << "...\n";
    m_input_mono_data = read_file_and_get_mono_data(infile, m_input_sfinfo, m_input_buffer);

    std::cout << "Loading reference track: " <<reffile << "...\n";
    m_ref_mono_data = read_file_and_get_mono_data(reffile, m_ref_sfinfo, m_ref_buffer);

    if (m_input_sfinfo.samplerate != m_ref_sfinfo.samplerate) {
        throw std::runtime_error("Error: Input and Reference files must have the same sample rate!");
    }
}

void FFTBalance::process_and_write(const std::string& outfile)
{
    m_bands = generate_bands(m_input_sfinfo.samplerate, m_num_bands_target);
    if (m_bands.empty()) {
        throw std::runtime_error("Error: Failed to generate valid frequency bands for the given sample rate.");
    }

    const int num_channels = m_input_sfinfo.channels;
    const sf_count_t num_frames = m_input_sfinfo.frames;
    const size_t N = num_frames;

    std::cout << "Spectral analysis parameters: Sample Rate=" << m_input_sfinfo.samplerate
              << " Hz, Channels=" << num_channels << ", Frames=" << num_frames
              << ", Bands=" << m_bands.size() << "\n";

    m_input_band_amp = calculate_band_amp(m_input_mono_data, m_input_sfinfo, m_bands);
    m_ref_band_amp = calculate_band_amp(m_ref_mono_data, m_ref_sfinfo, m_bands);

    std::cout << "Spectral analysis complete. Calculating gains...\n";
    calculate_gains();

    std::vector<std::vector<double>> processed_channels(num_channels, std::vector<double>(num_frames));
    double max_val_global = 0.0;

    m_channel_data.resize(num_frames);
    m_fft_channel = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1));
    m_ifft_channel = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    m_forward_channel_plan = fftw_plan_dft_r2c_1d(N, m_channel_data.data(), m_fft_channel, FFTW_ESTIMATE);
    m_backward_channel_plan = fftw_plan_dft_c2r_1d(N, m_ifft_channel, m_channel_data.data(), FFTW_ESTIMATE);

    std::cout << "Applying EQ to " << num_channels << " channels...\n";

    for (int c = 0; c < num_channels; ++c) {
        for (sf_count_t i = 0; i < num_frames; ++i) {
            m_channel_data[i] = m_input_buffer[i * num_channels + c];
        }

        fftw_execute(m_forward_channel_plan);

        for (size_t i = 0; i <= N / 2; i++) {
            m_ifft_channel[i][0] = m_fft_channel[i][0];
            m_ifft_channel[i][1] = m_fft_channel[i][1];
        }

        for (size_t b = 0; b < m_bands.size(); ++b) {
            double G = m_gains[b];
            size_t start = freqToBin(m_bands[b].low, m_input_sfinfo.samplerate, N);
            size_t end = freqToBin(m_bands[b].high, m_input_sfinfo.samplerate, N);

            if (end >= N / 2) end = N / 2;
            if (start >= end) continue;
            if (start == 0 && m_bands[b].low > 0) start = 1;

            for (size_t k = start; k <= end; ++k) {
                m_ifft_channel[k][0] *= G;
                m_ifft_channel[k][1] *= G;
            }
        }

        for (size_t k = 1; k < N / 2; ++k) {
            size_t nk = N - k;
            m_ifft_channel[nk][0] = m_ifft_channel[k][0];
            m_ifft_channel[nk][1] = -m_ifft_channel[k][1];
        }

        fftw_execute(m_backward_channel_plan);

        for (sf_count_t i = 0; i < num_frames; ++i) {
            double v = m_channel_data[i] / N;
            processed_channels[c][i] = v;
            if (std::abs(v) > max_val_global) max_val_global = std::abs(v);
        }
    }

    if (max_val_global < 1e-12) max_val_global = 1.0;

    const double final_scale = (max_val_global > 1.0) ? 1.0 / max_val_global : 1.0;
    std::cout << "Applying final output scaling factor: " << final_scale << "\n";

    std::vector<double> final_interleaved_data(num_frames * num_channels);
    for (sf_count_t i = 0; i < num_frames; ++i) {
        for (int c = 0; c < num_channels; ++c) {
            final_interleaved_data[i * num_channels + c] = processed_channels[c][i] * final_scale;
        }
    }

    SF_INFO outInfo = m_input_sfinfo;
    outInfo.channels = num_channels;

    SNDFILE* outFile = sf_open(outfile.c_str(), SFM_WRITE, &outInfo);
    if (!outFile) {
        throw std::runtime_error("Error opening output file: " + std::string(sf_strerror(NULL)));
    }

    sf_writef_double(outFile, final_interleaved_data.data(), num_frames);
    sf_close(outFile);

    std::cout << "Successfully matched spectral balance and wrote to: " << outfile << "\n";
}

std::vector<double> FFTBalance::read_file_and_get_mono_data(const std::string& filepath, SF_INFO& sfinfo_out, std::vector<double>& full_buffer_out)
{
    SNDFILE* inFile = sf_open(filepath.c_str(), SFM_READ, &sfinfo_out);
    if (!inFile) {
        throw std::runtime_error("Error opening file " + filepath + ": " + sf_strerror(NULL));
    }

    if (sfinfo_out.frames <= 0 || sfinfo_out.channels <= 0 || sfinfo_out.samplerate <= 0) {
        sf_close(inFile);
        throw std::runtime_error("Error: Invalid file info for " + filepath);
    }

    full_buffer_out.resize(sfinfo_out.frames * sfinfo_out.channels);
    sf_readf_double(inFile, full_buffer_out.data(), sfinfo_out.frames);
    sf_close(inFile);

    const sf_count_t N = sfinfo_out.frames;
    const int num_channels = sfinfo_out.channels;

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

std::vector<Band> FFTBalance::generate_bands(double samplerate, int num_bands)
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
        std::cerr << "Warning: Samplerate is too low to define full spectral range.\n";
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
            high = f_end;
        }
        double f_center = getBandCenterFrequency(low, high);
        bands.push_back({low, high, f_center});
        low = high;
    }
    return bands;
}

double FFTBalance::getBandCenterFrequency(double low, double high)
{
    return std::sqrt(low * high);
}

size_t FFTBalance::freqToBin(double freq, double fs, size_t N)
{
    if (N == 0 || fs == 0) return 0;
    return static_cast<size_t>(std::round(freq / fs * N));
}

std::vector<double> FFTBalance::calculate_band_amp(const std::vector<double>& mono_data, const SF_INFO& sfinfo_in, const std::vector<Band>& bands)
{
    const sf_count_t N = mono_data.size();
    if (N < 2) {
        throw std::runtime_error("Error: Sample size is too small for FFT.");
    }

    fftw_complex* fft_mono = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1));
    fftw_plan forward_mono = fftw_plan_dft_r2c_1d(N, const_cast<double*>(mono_data.data()), fft_mono, FFTW_ESTIMATE);
    fftw_execute(forward_mono);

    std::vector<double> band_amp;
    const double fs = sfinfo_in.samplerate;

    for (const auto& band : bands) {
        size_t start = freqToBin(band.low, fs, N);
        size_t end = freqToBin(band.high, fs, N);

        if (end >= N / 2) end = N / 2;
        if (start >= end) {
            band_amp.push_back(1e-12);
            continue;
        }

        if (start == 0 && band.low > 0) start = 1;

        double sum_sq_mag = 0;
        for (size_t k = start; k <= end; ++k) {
            double mag_sq = fft_mono[k][0] * fft_mono[k][0] + fft_mono[k][1] * fft_mono[k][1];
            sum_sq_mag += (k > 0 && k < N / 2) ? 2.0 * mag_sq : mag_sq;
        }
        band_amp.push_back(std::sqrt(sum_sq_mag));
    }

    fftw_destroy_plan(forward_mono);
    fftw_free(fft_mono);

    return band_amp;
}

void FFTBalance::calculate_gains()
{
    double input_overall_amp = std::accumulate(m_input_band_amp.begin(), m_input_band_amp.end(), 0.0) / m_input_band_amp.size();
    double ref_overall_amp = std::accumulate(m_ref_band_amp.begin(), m_ref_band_amp.end(), 0.0) / m_ref_band_amp.size();

    const double overall_max_amp = std::max(input_overall_amp, ref_overall_amp);
    const double AMPLITUDE_FLOOR = std::max(overall_max_amp * 1e-7, 1e-12);

    const double level_match_factor = (ref_overall_amp > 1e-12) ? input_overall_amp / ref_overall_amp : 1.0;

    std::cout << "Input Overall Avg Band RMS: " << input_overall_amp << "\n";
    std::cout << "Reference Overall Avg Band RMS: " << ref_overall_amp << "\n";
    std::cout << "Applying level match factor (Ref * " << level_match_factor << ") to normalize spectra before EQ.\n";
    std::cout << "Using Dynamic Amplitude Floor: " << AMPLITUDE_FLOOR << "\n";

    const double MAX_GAIN = std::pow(10.0, m_max_boost_db / 20.0);
    const double MIN_GAIN = std::pow(10.0, -m_max_cut_db / 20.0);

    std::cout << "Using Gain Limits: Boost=" << m_max_boost_db << " dB (factor " << MAX_GAIN
              << "), Cut=" << m_max_cut_db << " dB (factor " << MIN_GAIN << ").\n";

    for (size_t b = 0; b < m_bands.size(); ++b) {
        double target_amp_raw = m_ref_band_amp[b] * level_match_factor;
        double current_amp_floored = std::max(m_input_band_amp[b], AMPLITUDE_FLOOR);
        double target_amp_floored = std::max(target_amp_raw, AMPLITUDE_FLOOR);

        double G = target_amp_floored / current_amp_floored;

        if (G > MAX_GAIN) G = MAX_GAIN;
        if (G < MIN_GAIN) G = MIN_GAIN;

        m_gains.push_back(G);

        double dB = 20.0 * std::log10(G);
        std::cout << "Band " << b + 1 << " (" << m_bands[b].low << " Hz - " << m_bands[b].high
                  << " Hz): Gain = " << dB << " dB\n";
    }
}
