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
{
}

FFTBalance::~FFTBalance() = default;

void FFTBalance::loadTracks(const std::string & infile, const std::string & reffile)
{
    std::cout << "Loading input track: " << infile << "..." << std::endl;
    m_inputData = readFile(infile);

    std::cout << "Loading reference track: " << reffile << "..." << std::endl;
    m_referenceData = readFile(reffile);

    if (m_inputData.info.samplerate != m_referenceData.info.samplerate) {
        throw std::runtime_error("Error: Input and Reference files must have the same sample rate!");
    }
}

void FFTBalance::processAndWrite(const std::string & outfile)
{
    m_bands = generateBands(m_inputData.info.samplerate, m_bandCount);
    if (m_bands.empty()) {
        throw std::runtime_error("Error: Failed to generate valid frequency bands for the given sample rate.");
    }

    const size_t channelCount = m_inputData.info.channels;
    const size_t frameCount = m_inputData.info.frames;
    std::cout << "Spectral analysis parameters: Sample Rate=" << m_inputData.info.samplerate
              << " Hz, Channels=" << channelCount << ", Frames=" << frameCount
              << ", Bands=" << m_bands.size() << std::endl;

    // Determine the non-zero range for input data by checking all channels
    const double silenceThreshold = 1e-9;
    size_t inputFirstNonZeroFrame = frameCount;
    for (size_t i = 0; i < m_inputData.fullBuffer.size(); ++i) {
        if (std::abs(m_inputData.fullBuffer[i]) > silenceThreshold) {
            inputFirstNonZeroFrame = i / channelCount;
            break;
        }
    }

    size_t inputLastNonZeroFrame = 0;
    if (inputFirstNonZeroFrame < frameCount) { // If not all silent
        for (size_t i = m_inputData.fullBuffer.size(); i > 0; --i) {
            if (std::abs(m_inputData.fullBuffer[i - 1]) > silenceThreshold) {
                inputLastNonZeroFrame = (i - 1) / channelCount;
                break;
            }
        }
    }

    if (inputFirstNonZeroFrame >= frameCount) {
        std::cout << "Input track is completely silent. Skipping processing and writing silent output." << std::endl;
        SF_INFO outInfo { m_inputData.info };
        outInfo.channels = channelCount;
        if (const auto outFile = sf_open(outfile.c_str(), SFM_WRITE, &outInfo); !outFile) {
            throw std::runtime_error("Error opening output file: " + std::string(sf_strerror(nullptr)));
        } else {
            SampleVector silentData(frameCount * channelCount, 0.0);
            sf_writef_double(outFile, silentData.data(), frameCount);
            sf_close(outFile);
            std::cout << "Successfully wrote silent output to: " << outfile << std::endl;
        }
        return;
    }

    // Determine the non-zero range for reference mono data
    size_t refFirstNonZeroFrame = frameCount;
    size_t refLastNonZeroFrame = 0;
    for (size_t i = 0; i < frameCount; ++i) {
        if (m_referenceData.monoData[i] != 0.0) {
            if (refFirstNonZeroFrame == frameCount) {
                refFirstNonZeroFrame = i;
            }
            refLastNonZeroFrame = i;
        }
    }

    m_inputBandAmp = calculateBandAmp(m_inputData.monoData, m_inputData.info, m_bands, inputFirstNonZeroFrame, inputLastNonZeroFrame);
    m_refBandAmp = calculateBandAmp(m_referenceData.monoData, m_referenceData.info, m_bands, refFirstNonZeroFrame, refLastNonZeroFrame);

    std::cout << "Spectral analysis complete. Calculating gains..." << std::endl;
    calculateGains();

    std::vector<SampleVector> processedChannels(channelCount, SampleVector(frameCount, 0.0));
    double maxValGlobal { 0.0 };

    std::cout << "Applying EQ to " << channelCount << " channels...\n";

    for (size_t channelIndex = 0; channelIndex < channelCount; ++channelIndex) {
        // Copy silent part at the beginning from the original input
        for (size_t frameIndex = 0; frameIndex < inputFirstNonZeroFrame; ++frameIndex) {
            processedChannels[channelIndex][frameIndex] = m_inputData.fullBuffer[frameIndex * channelCount + channelIndex];
        }

        // Process the non-silent part using Weighted Overlap-Add (WOLA)
        const size_t nonSilentFrameCount = inputLastNonZeroFrame - inputFirstNonZeroFrame + 1;
        if (nonSilentFrameCount > 0) {
            const size_t chunkSize = 8192;
            const size_t hopSize = chunkSize / 2;
            SampleVector window(chunkSize);
            for (size_t i = 0; i < chunkSize; ++i) {
                window[i] = 0.5 * (1.0 - std::cos(2.0 * M_PI * i / (chunkSize - 1)));
            }

            SampleVector nonSilentData(nonSilentFrameCount);
            for (size_t i = 0; i < nonSilentFrameCount; ++i) {
                nonSilentData[i] = m_inputData.fullBuffer[(inputFirstNonZeroFrame + i) * channelCount + channelIndex];
            }

            size_t numChunks = (nonSilentFrameCount > 0) ? 1 + (nonSilentFrameCount - 1) / hopSize : 0;
            size_t paddedSize = (numChunks > 0) ? (numChunks - 1) * hopSize + chunkSize : 0;
            nonSilentData.resize(paddedSize, 0.0);

            SampleVector processedNonSilentData(paddedSize, 0.0);
            SampleVector windowSum(paddedSize, 0.0);

            size_t currentPos = 0;
            while (currentPos + chunkSize <= paddedSize) {
                SampleVector chunk(chunkSize);
                for (size_t i = 0; i < chunkSize; ++i) {
                    chunk[i] = nonSilentData[currentPos + i] * window[i];
                }

                auto * fftChannel = reinterpret_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * (chunkSize / 2 + 1)));
                auto forwardPlan = fftw_plan_dft_r2c_1d(chunkSize, chunk.data(), fftChannel, FFTW_ESTIMATE);
                fftw_execute(forwardPlan);

                for (size_t b = 0; b < m_bands.size(); ++b) {
                    const auto G = m_gains[b];
                    size_t start = freqToBin(m_bands[b].low, m_inputData.info.samplerate, chunkSize);
                    size_t end = freqToBin(m_bands[b].high, m_inputData.info.samplerate, chunkSize);
                    end = std::min(end, chunkSize / 2);
                    if (start < end) {
                        for (size_t k = start; k <= end; ++k) {
                            fftChannel[k][0] *= G;
                            fftChannel[k][1] *= G;
                        }
                    }
                }

                auto backwardPlan = fftw_plan_dft_c2r_1d(chunkSize, fftChannel, chunk.data(), FFTW_ESTIMATE);
                fftw_execute(backwardPlan);

                for (size_t i = 0; i < chunkSize; ++i) {
                    processedNonSilentData[currentPos + i] += (chunk[i] / chunkSize) * window[i];
                    windowSum[currentPos + i] += window[i] * window[i];
                }

                fftw_destroy_plan(forwardPlan);
                fftw_destroy_plan(backwardPlan);
                fftw_free(fftChannel);

                currentPos += hopSize;
            }

            for (size_t i = 0; i < nonSilentFrameCount; ++i) {
                if (windowSum[i] > 1e-9) {
                    processedNonSilentData[i] /= windowSum[i];
                }
                processedChannels[channelIndex][inputFirstNonZeroFrame + i] = processedNonSilentData[i];
            }
        }

        // Copy silent part at the end from the original input
        for (size_t frameIndex = inputLastNonZeroFrame + 1; frameIndex < frameCount; ++frameIndex) {
            processedChannels[channelIndex][frameIndex] = m_inputData.fullBuffer[frameIndex * channelCount + channelIndex];
        }

        for (size_t frameIndex = 0; frameIndex < frameCount; ++frameIndex) {
            if (std::abs(processedChannels[channelIndex][frameIndex]) > maxValGlobal) {
                maxValGlobal = std::abs(processedChannels[channelIndex][frameIndex]);
            }
        }
    }

    if (maxValGlobal < 1e-12) {
        maxValGlobal = 1.0;
    }

    const double finalScale = maxValGlobal > 1.0 ? 1.0 / maxValGlobal : 1.0;
    std::cout << "Applying final output scaling factor: " << finalScale << std::endl;

    SampleVector finalInterleavedData;
    finalInterleavedData.resize(frameCount * channelCount);
    for (size_t frameIndex = 0; frameIndex < frameCount; frameIndex++) {
        for (size_t channelIndex = 0; channelIndex < channelCount; channelIndex++) {
            finalInterleavedData[frameIndex * channelCount + channelIndex] = processedChannels[channelIndex][frameIndex] * finalScale;
        }
    }

    SF_INFO outInfo { m_inputData.info };
    outInfo.channels = channelCount;

    if (const auto outFile = sf_open(outfile.c_str(), SFM_WRITE, &outInfo); !outFile) {
        throw std::runtime_error("Error opening output file: " + std::string(sf_strerror(nullptr)));
    } else {
        sf_writef_double(outFile, finalInterleavedData.data(), frameCount);
        sf_close(outFile);
        std::cout << "Successfully matched spectral balance and wrote to: " << outfile << std::endl;
    }
}

FFTBalance::InputData FFTBalance::readFile(const std::string & filepath) const
{
    InputData inputData;
    if (const auto inFile = sf_open(filepath.c_str(), SFM_READ, &inputData.info); !inFile) {
        throw std::runtime_error("Error opening file " + filepath + ": " + sf_strerror(nullptr));
    } else {
        if (inputData.info.frames <= 0 || inputData.info.channels <= 0 || inputData.info.samplerate <= 0) {
            sf_close(inFile);
            throw std::runtime_error("Error: Invalid file info for " + filepath);
        }
        inputData.fullBuffer.resize(inputData.info.frames * inputData.info.channels);
        sf_readf_double(inFile, inputData.fullBuffer.data(), inputData.info.frames);
        sf_close(inFile);

        const size_t frameCount = inputData.info.frames;
        const size_t channelCount = inputData.info.channels;
        inputData.monoData.resize(frameCount);
        for (size_t frameIndex = 0; frameIndex < frameCount; frameIndex++) {
            double sum = 0;
            for (size_t channelIndex = 0; channelIndex < channelCount; channelIndex++) {
                sum += inputData.fullBuffer.at(frameIndex * channelCount + channelIndex);
            }
            inputData.monoData.at(frameIndex) = sum / static_cast<double>(channelCount);
        }
        return inputData;
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
        std::cerr << "Warning: Samplerate is too low to define full spectral range." << std::endl;
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

        // Rounding logic for nice numbers
        if (low < 100.0) {
            low = std::round(low);
        } else if (low < 1000.0) {
            low = std::round(low / 10.0) * 10.0;
        } else {
            low = std::round(low / 100.0) * 100.0;
        }

        if (high < 100.0) {
            high = std::round(high);
        } else if (high < 1000.0) {
            high = std::round(high / 10.0) * 10.0;
        } else {
            high = std::round(high / 100.0) * 100.0;
        }

        // Ensure high is always greater than low after rounding
        if (high <= low) {
            high = low + 1.0; // Smallest possible increment
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

FFTBalance::SampleVector FFTBalance::calculateBandAmp(
    const SampleVector & monoData,
    const SF_INFO & sfinfoIn,
    const std::vector<Band> & bands,
    size_t firstNonZeroFrame,
    size_t lastNonZeroFrame) const
{
    if (firstNonZeroFrame > lastNonZeroFrame || monoData.empty()) {
        return SampleVector(bands.size(), 1e-12);
    }

    const size_t nonSilentFrameCount = lastNonZeroFrame - firstNonZeroFrame + 1;
    if (nonSilentFrameCount < 2) {
        return SampleVector(bands.size(), 1e-12);
    }

    const size_t chunkSize = 8192;
    const size_t hopSize = chunkSize / 2;
    SampleVector window(chunkSize);
    for (size_t i = 0; i < chunkSize; ++i) {
        window[i] = 0.5 * (1.0 - std::cos(2.0 * M_PI * i / (chunkSize - 1)));
    }

    SampleVector nonSilentData(nonSilentFrameCount);
    for (size_t i = 0; i < nonSilentFrameCount; ++i) {
        nonSilentData[i] = monoData[firstNonZeroFrame + i];
    }

    SampleVector bandAmpsSum(bands.size(), 0.0);

    size_t currentPos = 0;
    size_t chunksProcessed = 0;
    while (currentPos + chunkSize <= nonSilentFrameCount) {
        SampleVector chunk(chunkSize);
        for (size_t i = 0; i < chunkSize; ++i) {
            chunk[i] = nonSilentData[currentPos + i] * window[i];
        }

        auto * fftMono = reinterpret_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * (chunkSize / 2 + 1)));
        auto forwardMono = fftw_plan_dft_r2c_1d(chunkSize, chunk.data(), fftMono, FFTW_ESTIMATE);
        fftw_execute(forwardMono);

        for (size_t b = 0; b < bands.size(); ++b) {
            const auto & band = bands[b];
            auto startBin = freqToBin(band.low, sfinfoIn.samplerate, chunkSize);
            auto endBin = freqToBin(band.high, sfinfoIn.samplerate, chunkSize);
            endBin = std::min(endBin, chunkSize / 2);

            double sumSqMag = 0.0;
            if (startBin < endBin) {
                for (size_t k = startBin; k <= endBin; ++k) {
                    const double magSq = fftMono[k][0] * fftMono[k][0] + fftMono[k][1] * fftMono[k][1];
                    sumSqMag += magSq;
                }
            }
            bandAmpsSum[b] += std::sqrt(sumSqMag);
        }

        fftw_destroy_plan(forwardMono);
        fftw_free(fftMono);
        currentPos += hopSize;
        chunksProcessed++;
    }

    if (chunksProcessed == 0) {
        return SampleVector(bands.size(), 1e-12);
    }

    SampleVector avgBandAmp(bands.size());
    for (size_t b = 0; b < bands.size(); ++b) {
        avgBandAmp[b] = bandAmpsSum[b] / chunksProcessed;
    }

    return avgBandAmp;
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
