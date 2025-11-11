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

#include "fft_balance.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>

FFTBalance::FFTBalance(int bandCount, double maxBoostDb, double maxCutDb)
  : m_bandCount { bandCount }
  , m_maxBoostDb { maxBoostDb }
  , m_maxCutDb { maxCutDb }
  , m_forwardChannelPlan { nullptr }
  , m_backwardChannelPlan { nullptr }
  , m_fftChannel { nullptr }
  , m_ifftChannel { nullptr }
{
}

FFTBalance::~FFTBalance()
{
    if (m_forwardChannelPlan) {
        fftw_destroy_plan(m_forwardChannelPlan);
    }
    if (m_backwardChannelPlan) {
        fftw_destroy_plan(m_backwardChannelPlan);
    }
    if (m_fftChannel) {
        fftw_free(m_fftChannel);
    }
    if (m_ifftChannel) {
        fftw_free(m_ifftChannel);
    }
}

void FFTBalance::loadTracks(const std::string & infile, const std::string & reffile)
{
    std::cout << "Loading input track: " << infile << "...\n";
    m_inputMonoData = readFileAndGetMonoData(infile, m_inputSfinfo, m_inputBuffer);

    std::cout << "Loading reference track: " << reffile << "...\n";
    m_refMonoData = readFileAndGetMonoData(reffile, m_refSfinfo, m_refBuffer);

    if (m_inputSfinfo.samplerate != m_refSfinfo.samplerate) {
        throw std::runtime_error("Error: Input and Reference files must have the same sample rate!");
    }
}

void FFTBalance::processAndWrite(const std::string & outfile)
{
    m_bands = generateBands(m_inputSfinfo.samplerate, m_bandCount);
    if (m_bands.empty()) {
        throw std::runtime_error("Error: Failed to generate valid frequency bands for the given sample rate.");
    }

    const size_t channelCount = m_inputSfinfo.channels;
    const size_t frameCount = m_inputSfinfo.frames;
    std::cout << "Spectral analysis parameters: Sample Rate=" << m_inputSfinfo.samplerate
              << " Hz, Channels=" << channelCount << ", Frames=" << frameCount
              << ", Bands=" << m_bands.size() << std::endl;

    m_inputBandAmp = calculateBandAmp(m_inputMonoData, m_inputSfinfo, m_bands);
    m_refBandAmp = calculateBandAmp(m_refMonoData, m_refSfinfo, m_bands);

    std::cout << "Spectral analysis complete. Calculating gains..." << std::endl;
    calculateGains();

    std::vector<std::vector<double>> processedChannels(channelCount, std::vector<double>(frameCount));
    double maxValGlobal { 0.0 };

    m_channelData.resize(frameCount);
    m_fftChannel = reinterpret_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * (frameCount / 2 + 1)));
    m_ifftChannel = reinterpret_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * frameCount));
    m_forwardChannelPlan = fftw_plan_dft_r2c_1d(frameCount, m_channelData.data(), m_fftChannel, FFTW_ESTIMATE);
    m_backwardChannelPlan = fftw_plan_dft_c2r_1d(frameCount, m_ifftChannel, m_channelData.data(), FFTW_ESTIMATE);

    std::cout << "Applying EQ to " << channelCount << " channels...\n";

    for (size_t channelIndex = 0; channelIndex < channelCount; ++channelIndex) {
        for (size_t frameIndex = 0; frameIndex < frameCount; ++frameIndex) {
            m_channelData[frameIndex] = m_inputBuffer[frameIndex * channelCount + channelIndex];
        }

        fftw_execute(m_forwardChannelPlan);

        for (size_t i = 0; i <= frameCount / 2; i++) {
            m_ifftChannel[i][0] = m_fftChannel[i][0];
            m_ifftChannel[i][1] = m_fftChannel[i][1];
        }

        for (size_t b = 0; b < m_bands.size(); ++b) {
            const auto G = m_gains[b];
            size_t start = freqToBin(m_bands[b].low, m_inputSfinfo.samplerate, frameCount);
            size_t end = freqToBin(m_bands[b].high, m_inputSfinfo.samplerate, frameCount);
            end = std::min(end, frameCount / 2);
            if (start < end) {
                if (start == 0 && m_bands[b].low > 0) {
                    start = 1;
                }
                for (size_t k = start; k <= end; ++k) {
                    m_ifftChannel[k][0] *= G;
                    m_ifftChannel[k][1] *= G;
                }
            }
        }

        for (size_t frameIndex = 1; frameIndex < frameCount / 2; ++frameIndex) {
            size_t mirroredIndex = frameCount - frameIndex;
            m_ifftChannel[mirroredIndex][0] = m_ifftChannel[frameIndex][0];
            m_ifftChannel[mirroredIndex][1] = -m_ifftChannel[frameIndex][1];
        }

        fftw_execute(m_backwardChannelPlan);

        for (size_t frameIndex = 0; frameIndex < frameCount; ++frameIndex) {
            double v = m_channelData[frameIndex] / frameCount;
            processedChannels[channelIndex][frameIndex] = v;
            if (std::abs(v) > maxValGlobal)
                maxValGlobal = std::abs(v);
        }
    }

    if (maxValGlobal < 1e-12) {
        maxValGlobal = 1.0;
    }

    const double finalScale = maxValGlobal > 1.0 ? 1.0 / maxValGlobal : 1.0;
    std::cout << "Applying final output scaling factor: " << finalScale << std::endl;

    std::vector<double> finalInterleavedData;
    finalInterleavedData.resize(frameCount * channelCount);
    for (size_t frameIndex = 0; frameIndex < frameCount; frameIndex++) {
        for (size_t channelIndex = 0; channelIndex < channelCount; channelIndex++) {
            finalInterleavedData[frameIndex * channelCount + channelIndex] = processedChannels[channelIndex][frameIndex] * finalScale;
        }
    }

    SF_INFO outInfo { m_inputSfinfo };
    outInfo.channels = channelCount;

    if (const auto outFile = sf_open(outfile.c_str(), SFM_WRITE, &outInfo); !outFile) {
        throw std::runtime_error("Error opening output file: " + std::string(sf_strerror(nullptr)));
    } else {
        sf_writef_double(outFile, finalInterleavedData.data(), frameCount);
        sf_close(outFile);
        std::cout << "Successfully matched spectral balance and wrote to: " << outfile << std::endl;
    }
}

std::vector<double> FFTBalance::readFileAndGetMonoData(const std::string & filepath, SF_INFO & sfinfo_out, std::vector<double> & full_buffer_out)
{
    if (const auto inFile = sf_open(filepath.c_str(), SFM_READ, &sfinfo_out); !inFile) {
        throw std::runtime_error("Error opening file " + filepath + ": " + sf_strerror(nullptr));
    } else {
        if (sfinfo_out.frames <= 0 || sfinfo_out.channels <= 0 || sfinfo_out.samplerate <= 0) {
            sf_close(inFile);
            throw std::runtime_error("Error: Invalid file info for " + filepath);
        }
        full_buffer_out.resize(sfinfo_out.frames * sfinfo_out.channels);
        sf_readf_double(inFile, full_buffer_out.data(), sfinfo_out.frames);
        sf_close(inFile);

        const size_t frameCount = sfinfo_out.frames;
        const size_t channelCount = sfinfo_out.channels;
        std::vector<double> monoData;
        monoData.resize(frameCount);
        for (size_t frameIndex = 0; frameIndex < frameCount; frameIndex++) {
            double sum = 0;
            for (size_t channelIndex = 0; channelIndex < channelCount; channelIndex++) {
                sum += full_buffer_out[frameIndex * channelCount + channelIndex];
            }
            monoData.at(frameIndex) = sum / static_cast<double>(channelCount);
        }
        return monoData;
    }
}

std::vector<Band> FFTBalance::generateBands(double samplerate, size_t bandCount) const
{
    if (bandCount < 1) {
        throw std::runtime_error { "Error: Number of bands must be at least 1." };
    }

    const double startFreq = 20.0;
    const double nyquist = samplerate / 2.0;
    const double endFreq = std::min(20000.0, nyquist - 1.0);

    std::vector<Band> bands;
    if (endFreq <= startFreq) {
        std::cerr << "Warning: Samplerate is too low to define full spectral range.\n";
        if (nyquist > 20.0) {
            bands.push_back({ startFreq, nyquist, getBandCenterFrequency(startFreq, nyquist) });
        }
        return bands;
    }

    const double multiplier = std::pow(endFreq / startFreq, 1.0 / static_cast<int>(bandCount));
    double low = startFreq;
    for (size_t bandIndex = 0; bandIndex < bandCount; bandIndex++) {
        double high = low * multiplier;
        if (bandIndex == bandCount - 1) {
            high = endFreq;
        }
        const double centerFreq { getBandCenterFrequency(low, high) };
        bands.push_back({ low, high, centerFreq });
        low = high;
    }
    return bands;
}

double FFTBalance::getBandCenterFrequency(double low, double high) const
{
    return std::sqrt(low * high);
}

size_t FFTBalance::freqToBin(double frequency, double samplerate, size_t binCount) const
{
    if (binCount == 0 || samplerate < 0.001) {
        return 0;
    } else {
        return static_cast<size_t>(std::round(frequency / samplerate * static_cast<int>(binCount)));
    }
}

std::vector<double> FFTBalance::calculateBandAmp(const std::vector<double> & monoData, const SF_INFO & sfinfoIn, const std::vector<Band> & bands) const
{
    const auto N = monoData.size();
    if (N < 2) {
        throw std::runtime_error { "Error: Sample size is too small for FFT." };
    }

    const auto fftMono = reinterpret_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1)));
    const auto forwardMono = fftw_plan_dft_r2c_1d(N, const_cast<double *>(monoData.data()), fftMono, FFTW_ESTIMATE);
    fftw_execute(forwardMono);

    std::vector<double> bandAmp;
    const auto samplerate = sfinfoIn.samplerate;
    for (const auto & band : bands) {
        auto startBin = freqToBin(band.low, samplerate, N);
        auto endBin = freqToBin(band.high, samplerate, N);
        endBin = std::min(endBin, N / 2);
        if (startBin >= endBin) {
            bandAmp.push_back(1e-12);
        } else {
            if (startBin == 0 && band.low > 0) {
                startBin = 1;
            }
            double sumSqMag { 0.0 };
            for (size_t k = startBin; k <= endBin; ++k) {
                const double magSq = fftMono[k][0] * fftMono[k][0] + fftMono[k][1] * fftMono[k][1];
                sumSqMag += (k > 0 && k < N / 2) ? 2.0 * magSq : magSq;
            }
            bandAmp.push_back(std::sqrt(sumSqMag));
        }
    }

    fftw_destroy_plan(forwardMono);
    fftw_free(fftMono);

    return bandAmp;
}

void FFTBalance::calculateGains()
{
    const double inputOverallAmp = std::accumulate(m_inputBandAmp.begin(), m_inputBandAmp.end(), 0.0) / m_inputBandAmp.size();
    const double refOverallAmp = std::accumulate(m_refBandAmp.begin(), m_refBandAmp.end(), 0.0) / m_refBandAmp.size();

    const double overallMaxAmp = std::max(inputOverallAmp, refOverallAmp);
    const double amplitudeFloor = std::max(overallMaxAmp * 1e-7, 1e-12);

    const double levelMatchFactor = refOverallAmp > 1e-12 ? inputOverallAmp / refOverallAmp : 1.0;

    std::cout << "Input Overall Avg Band RMS: " << inputOverallAmp << std::endl;
    std::cout << "Reference Overall Avg Band RMS: " << refOverallAmp << std::endl;
    std::cout << "Applying level match factor (Ref * " << levelMatchFactor << ") to normalize spectra before EQ." << std::endl;
    std::cout << "Using Dynamic Amplitude Floor: " << amplitudeFloor << std::endl;

    const double maxGain = std::pow(10.0, m_maxBoostDb / 20.0);
    const double minGain = std::pow(10.0, -m_maxCutDb / 20.0);

    std::cout << "Using Gain Limits: Boost=" << m_maxBoostDb << " dB (factor " << maxGain
              << "), Cut=" << m_maxCutDb << " dB (factor " << minGain << ")." << std::endl;

    for (size_t b = 0; b < m_bands.size(); ++b) {
        const double targetAmpRaw = m_refBandAmp[b] * levelMatchFactor;
        const double currentAmpFloored = std::max(m_inputBandAmp[b], amplitudeFloor);
        const double targetAmpFloored = std::max(targetAmpRaw, amplitudeFloor);

        double G = targetAmpFloored / currentAmpFloored;
        G = std::min(G, maxGain);
        G = std::max(G, minGain);
        m_gains.push_back(G);

        const double dB = 20.0 * std::log10(G);
        std::cout << "Band " << b + 1 << " (" << m_bands[b].low << " Hz - " << m_bands[b].high
                  << " Hz): Gain = " << dB << " dB" << std::endl;
    }
}
