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

#include "application.hpp"
#include "fft_balance.hpp"
#include <iostream>
#include <stdexcept>
#include <vector>

Application::Application(int argc, char ** argv)
  : m_numBands { 8 }
  , m_maxBoostDb { 3.0 }
  , m_maxCutDb { 3.0 }
{
    const std::vector<std::string> args { argv, argv + argc };
    if (argc < 4 || argc > 7) {
        printUsage(args.at(0));
        throw std::runtime_error("Invalid number of arguments.");
    }

    try {
        m_inFile = args.at(1);
        m_refFile = args.at(2);
        m_outFile = args.at(3);

        if (argc >= 5) {
            m_numBands = std::stoi(args.at(4));
            if (m_numBands < 2 || m_numBands > 30) {
                std::cerr << "Warning: Number of bands should be between 2 and 30. Using default 8." << std::endl;
                m_numBands = 8;
            }
        }

        if (argc >= 6) {
            m_maxBoostDb = std::stod(args.at(5));
            if (m_maxBoostDb < 0.0) {
                std::cerr << "Warning: Max boost dB cannot be negative. Using default 3.0." << std::endl;
                m_maxBoostDb = 3.0;
            }
        }

        if (argc == 7) {
            m_maxCutDb = std::stod(args.at(6));
            if (m_maxCutDb < 0.0) {
                std::cerr << "Warning: Max cut dB cannot be negative. Using default 3.0." << std::endl;
                m_maxCutDb = 3.0;
            }
        }
    } catch (const std::out_of_range &) {
        printUsage(args.at(0));
        throw std::runtime_error("Missing arguments.");
    } catch (const std::exception &) {
        throw std::runtime_error("Invalid numeric argument provided for optional parameters.");
    }
}

void Application::printUsage(const std::string applicationName) const
{
    std::cout << "Usage: " << applicationName << " input.wav reference.wav output.wav [num_bands (default 8)] "
                                                 "[max_boost_db (default 3.0)] [max_cut_db (default 3.0)]"
              << std::endl;
    std::cout << std::endl
              << "This tool balances the N-band spectral RMS of input.wav to match reference.wav." << std::endl;
}

int Application::run()
{
    FFTBalance balancer { m_numBands, m_maxBoostDb, m_maxCutDb };
    balancer.loadTracks(m_inFile, m_refFile);
    balancer.processAndWrite(m_outFile);
    return EXIT_SUCCESS;
}
