#include "scale.hpp"


extern "C" {

void scale_kernel(
    real_t * inout,
    const real_t mult,
    int size)
{
#pragma HLS INTERFACE m_axi port = inout bundle = gmem0 depth = 4096
    for (int i = 0; i < size; i++)
    {
        inout[i] *= mult;
    }
}
}
