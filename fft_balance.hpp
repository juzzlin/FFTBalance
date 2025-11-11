#ifndef FFTBALANCE_HPP
#define FFTBALANCE_HPP

#include <string>
#include <vector>
#include <sndfile.h>
#include <fftw3.h>

struct Band
{
    double low;
    double high;
    double f_center;
};

class FFTBalance
{
public:
    FFTBalance(int num_bands, double max_boost_db, double max_cut_db);
    ~FFTBalance();

    void load_tracks(const std::string& infile, const std::string& reffile);
    void process_and_write(const std::string& outfile);

private:
    // Helper functions
    std::vector<double> read_file_and_get_mono_data(const std::string& filepath, SF_INFO& sfinfo, std::vector<double>& full_buffer);
    std::vector<Band> generate_bands(double samplerate, int num_bands);
    double getBandCenterFrequency(double low, double high);
    size_t freqToBin(double freq, double fs, size_t N);
    std::vector<double> calculate_band_amp(const std::vector<double>& mono_data, const SF_INFO& sfinfo, const std::vector<Band>& bands);
    void calculate_gains();

    // Member variables
    int m_num_bands_target;
    double m_max_boost_db;
    double m_max_cut_db;

    SF_INFO m_input_sfinfo;
    SF_INFO m_ref_sfinfo;

    std::vector<double> m_input_buffer;
    std::vector<double> m_ref_buffer;
    std::vector<double> m_input_mono_data;
    std::vector<double> m_ref_mono_data;

    std::vector<Band> m_bands;
    std::vector<double> m_input_band_amp;
    std::vector<double> m_ref_band_amp;
    std::vector<double> m_gains;

    fftw_plan m_forward_channel_plan;
    fftw_plan m_backward_channel_plan;
    fftw_complex* m_fft_channel;
    fftw_complex* m_ifft_channel;
    std::vector<double> m_channel_data;
};

#endif // FFTBALANCE_HPP
