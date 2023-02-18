#ifndef HOST_RANDOM_H
#define HOST_RANDOM_H

#include <random>

/// Uniform random number generator
class rand_t {
    std::mt19937 gen_;
    std::uniform_real_distribution<> dist_;


public:
    template<typename T>
    rand_t(real_t a, real_t b, T seed = static_cast<real_t>(0)) :
      gen_(seed), dist_(-0.5, 0.5)
    {}

    real_t operator()()
    { return dist_(gen_); }
};

#endif
