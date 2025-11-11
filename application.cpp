#include "application.hpp"
#include "fft_balance.hpp"
#include <iostream>
#include <stdexcept>

Application::Application(int argc, char** argv)
    : m_num_bands(8), m_max_boost_db(3.0), m_max_cut_db(3.0), m_valid_args(false)
{
    if (argc < 4 || argc > 7) {
        print_usage();
        return;
    }

    m_infile = argv[1];
    m_reffile = argv[2];
    m_outfile = argv[3];

    try {
        if (argc >= 5) {
            m_num_bands = std::stoi(argv[4]);
            if (m_num_bands < 2 || m_num_bands > 30) {
                std::cerr << "Warning: Number of bands should be between 2 and 30. Using default 8.\n";
                m_num_bands = 8;
            }
        }

        if (argc >= 6) {
            m_max_boost_db = std::stod(argv[5]);
            if (m_max_boost_db < 0.0) {
                std::cerr << "Warning: Max boost dB cannot be negative. Using default 3.0.\n";
                m_max_boost_db = 3.0;
            }
        }

        if (argc == 7) {
            m_max_cut_db = std::stod(argv[6]);
            if (m_max_cut_db < 0.0) {
                std::cerr << "Warning: Max cut dB cannot be negative. Using default 3.0.\n";
                m_max_cut_db = 3.0;
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Warning: Invalid numeric argument provided for optional parameters. Using defaults.\n";
        m_num_bands = 8;
        m_max_boost_db = 3.0;
        m_max_cut_db = 3.0;
    }
    m_valid_args = true;
}

void Application::print_usage()
{
    std::cout << "Usage: fft_balance input.wav reference.wav output.wav [num_bands (default 8)] "
                 "[max_boost_db (default 3.0)] [max_cut_db (default 3.0)]\n";
    std::cout << "This tool balances the N-band spectral RMS of input.wav to match reference.wav.\n";
}

int Application::run()
{
    if (!m_valid_args) {
        return 1;
    }

    try {
        FFTBalance balancer(m_num_bands, m_max_boost_db, m_max_cut_db);
        balancer.load_tracks(m_infile, m_reffile);
        balancer.process_and_write(m_outfile);
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << "\n";
        return 1;
    } catch (...) {
        std::cerr << "An unknown error occurred.\n";
        return 1;
    }

    return 0;
}
