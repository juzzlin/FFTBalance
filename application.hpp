#ifndef APPLICATION_HPP
#define APPLICATION_HPP

#include <string>

class Application
{
public:
    Application(int argc, char** argv);
    int run();

private:
    void print_usage();

    std::string m_infile;
    std::string m_reffile;
    std::string m_outfile;
    int m_num_bands;
    double m_max_boost_db;
    double m_max_cut_db;
    bool m_valid_args;
};

#endif // APPLICATION_HPP
