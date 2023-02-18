#include "input.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <string>
#include <vector>

// trim from start (in place)
void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
void trim(std::string &s) {
    rtrim(s);
    ltrim(s);
}

std::vector<std::string> split(const std::string & str)
{   
    std::vector<std::string> toks;
 
    std::istringstream instr(str);
    std::string tok;
    while (std::getline(instr, tok, ' '))
    {
        if (!tok.empty()) toks.emplace_back(tok);
    }
    return toks;
}


input_t::input_t() : 
    dims{2, 2, 2},
    unit{1, 1, 1}
{}

input_t::input_t(const char * filename)
{
    if (!filename) return;

    std::ifstream infile(filename);
    std::string line;

    int lineno = 0;
    while (std::getline(infile, line))
    {
        ++lineno;

        auto found = line.find_first_of("#");
        if (found != std::string::npos) line.erase(found, line.size());
        
        auto splited = split(line);

        if (line.size() > 1) {
            auto ret = set_value(splited[0], splited.begin()+1, splited.end());
            if (ret) {
                std::cout << "Error on line " << lineno << " of " << filename << std::endl;
                exit(1);
            }            
        }

    }
 
}
            
std::ostream& operator<<(std::ostream & out, const input_t & o)
{
    out << "Temperature = " << o.temperature << std::endl;
    out << "Boltzmann   = " << o.boltzmann << std::endl;
    out << "Mass        = " << o.mass << std::endl;
    out << "dims        = " << o.dims[0] << ", " << o.dims[1] << ", " << o.dims[2]
        << std::endl;
    out << "unit length = " << o.unit[0] << ", " << o.unit[1] << ", "
        << o.unit[2] << std::endl;
    out << "Time step   = " << o.dt << std::endl;
    out << "Start time  = " << o.start_time << std::endl;
    out << "Num steps   = " << o.num_steps << std::endl;
    out << "Output freq = " << o.output_freq << std::endl;
    out << "Stats freq = " << o.stats_freq;
    return out;
}
