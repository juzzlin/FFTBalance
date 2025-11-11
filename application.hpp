#ifndef APPLICATION_HPP
#define APPLICATION_HPP

#include <string>

class Application
{
public:
    Application(int argc, char ** argv);

    int run();

private:
    void printUsage(const std::string applicationName) const;

    std::string m_inFile;
    std::string m_refFile;
    std::string m_outFile;

    int m_numBands;

    double m_maxBoostDb;
    double m_maxCutDb;
};

#endif // APPLICATION_HPP
