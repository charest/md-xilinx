#ifndef KERNEL_MOVE_H
#define KERNEL_MOVE_H

#include "../../config.h"
#include "hls_vector.h"

extern "C" {
void move_kernel(
    real_t * xinout,
    real_t * vin,
    real_t * ain,
    const real_t dt,
    const real_t box_length[3],
    const real_t mass,
    const real_t temperature,
    const real_t boltzmann,
    int size);
}

#endif
