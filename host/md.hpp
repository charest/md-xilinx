#ifndef HOST_MD_H
#define HOST_MD_H

#include "config.hpp"

class input_t;

class md_t {

    int dims_[3];
    real_t unit_[3];
    real_t box_length_[3];

    real_t boltzmann_;
    real_t mass_;
    real_t temperature_;

public:
    
    md_t(const input_t & in);

    int size_lattice() const;

    real_t init_particles(
        real_t * x,
        real_t * v,
        real_t * a,
        int num_parts) const;
    
    real_t move_particles(
        real_t * x,
        real_t * v,
        real_t * a,
        real_t dt,
        int num_parts) const;

};

#endif
