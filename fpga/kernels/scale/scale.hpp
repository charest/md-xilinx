#ifndef KERNEL_SCALE_H
#define KERNEL_SCALE_H

#include "../../config.h"

extern "C" {
void scale_kernel(
    real_t * inout,
    const real_t mult,
    int size);
}

#endif
