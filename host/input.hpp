#ifndef HOST_INPUT_H
#define HOST_INPUT_H

#include "config.hpp"

#include <iostream>
#include <string>

struct input_t {

    // temperature
    real_t temperature = 10;
    // boltzmann constant
    real_t boltzmann = 1;
    // mass
    real_t mass = 1;
    // grid
    int dims[3];
    // box size
    real_t unit[3];
    // time step
    real_t dt = 1.e-3;
    // number of steps
    int num_steps = 10;
    // start time
    real_t start_time = 0;
    // output frequency
    int output_freq = 0;
    int stats_freq = 1;

    input_t();
    input_t(const char *);

    template<typename It>
    int set_value(const std::string & key, It start, It end);
    
    template<typename T>
    int set_value(const char * str, T & var, const std::string & val);
    
    template<typename T, typename It>
    int set_value(const char * str, T * vbeg, T * vend, It beg, It end);
    
    friend std::ostream& operator<<(std::ostream&, const input_t & o);

};
    
template<>
inline int input_t::set_value<int>(const char * str, int & var, const std::string & val)
{
    var = std::stoi(val);
    return 0;
}

template<>
inline int input_t::set_value<real_t>(const char * str, real_t & var, const std::string & val)
{
    var = std::stod(val);
    return 0;
}
    
template<typename T, typename It>
inline int input_t::set_value(const char * str, T * vbeg, T * vend, It beg, It end)
{
    auto vlen = std::distance(vbeg, vend);
    auto ilen = std::distance(beg, end);
    if (ilen < vlen) {
        std::cout << "Variable '" << str << "' requires " << vlen << " values";
        std::cout << ", only " << ilen << " provided." << std::endl;
        return 1;
    }
    for (int i=0; i<vlen; ++i) {
        auto ret = set_value(str, vbeg[i], *(beg+i));
        if (ret) return ret;
    }
    return 0;
}

template<typename It>
inline int input_t::set_value(const std::string & key, It start, It end)
{           
    int ret = 0; 
    if (key == "temperature") {
        ret = set_value<real_t>("temperature", temperature, *start);
    }
    else if (key == "boltzmann") {
        ret = set_value<real_t>("boltzmann", boltzmann, *start);
    }
    else if (key == "mass") {
        ret = set_value<real_t>("mass", mass, *start);
    }
    else if (key == "dim") {
        ret = set_value<int>("dim", dims, dims+3, start, end);
    }
    else if (key == "box") {
        ret = set_value<real_t>("unit", unit, unit+3, start, end);
    }
    else if (key == "dt") {
        ret = set_value<real_t>("dt", dt, *start);
    }
    else if (key == "num_steps") {
        ret = set_value<int>("num_steps", num_steps, *start);
    }
    else if (key == "start_time") {
        ret = set_value<real_t>("start_time", start_time, *start);
    }
    else if (key == "output_freq") {
        ret = set_value<int>("output_freq", output_freq, *start);
    }
    else if (key == "stats_freq") {
        ret = set_value<int>("stats_freq", stats_freq, *start);
    }
    else {
        ret = 1;
        std::cout << "Unknown variable '" << key << "'" << std::endl;
    }
    return ret;
}


#endif
