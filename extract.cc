#include <fstream>
#include <iostream>
#include <string>
#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>

// trim from start (in place)
static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

int main(int argc, char *argv[])
{
    if (argc < 2)
        return -1;

    std::string argv0 = std::string(argv[1]);
    size_t startPos = std::stol(argv0);
    if (startPos < 1)
        return -1;
    
    std::cerr << "startPos:" << startPos << std::endl;

    const char *filename = "data/chromatid20.fa";
    size_t length = 150;

    std::string whole;

    std::ifstream ifs(filename);
    std::string str;
    if (ifs.fail())
    {
        std::cerr << "失敗" << std::endl;
        return -1;
    }

    int lines = 0;

    getline(ifs, str); // skip header
    while (getline(ifs, str))
    {
        if (++lines % 10000 == 0)
            std::cerr << ".";

        trim(str);
        whole += str;
    }
    std::cerr << std::endl;

    std::cout << whole.substr(startPos - 1, length) << std::endl;

    return 0;
}
