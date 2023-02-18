#include "move.hpp"

#include <cmath>

extern "C" {

void move_kernel(
    real_t * x,
    real_t * v,
    real_t * a,
    const real_t dt,
    const real_t box_length[3],
    const real_t mass,
    const real_t temperature,
    const real_t boltzmann,
    int size)
{
#pragma HLS INTERFACE m_axi port = v bundle = gmem0 depth = 4096
#pragma HLS INTERFACE m_axi port = a bundle = gmem1 depth = 4096
#pragma HLS INTERFACE m_axi port = x bundle = gmem2 depth = 4096

    //--- move particles

    auto dt_sqr = dt*dt;

    for (int i = 0; i < size; ++i)
    {
        for (int d=0; d<3; ++d)
        {
            // move to new position
            auto pos = 3*i + d;
            auto xd = x[pos] + v[pos]*dt + 0.5*dt_sqr*a[pos];
        
            // apply boundary condition
            auto len = box_length[d];
            if (xd > len)
                x[pos] = xd - len;
            else if (xd < 0)
                x[pos] = xd + len;
            else
                x[pos] = xd;
        }

    }

    //--- compute forces and update velocity/acceleration

    real_t en = 0;
    real_t vsum[] = {0, 0, 0};
    auto inv_mass = static_cast<real_t>(1) / mass;
    
    for (int i=0; i<size; ++i) {

        //--- sum all forces for particle i
        real_t f[] = {0, 0, 0};

        for (int j=0; j<size; ++j)
        {
            if (i != j) {
                // distance
                auto posi = 3*i;
                auto posj = 3*j;
                real_t xr[3] = {
                    x[posi + 0] - x[posj + 0],
                    x[posi + 1] - x[posj + 1],
                    x[posi + 2] - x[posj + 2],
                };
                // minimum image criterion
                for (int d=0; d<3; ++d)
                    xr[d] -= std::round(xr[d] / box_length[d]) * box_length[d];
                // dot product
                real_t r2 = 0;
                for (int d=0; d<3; ++d)
                    r2 += xr[d] * xr[d];
                r2 = 1 / r2;
                auto r6 = std::pow(r2, 3);
                // force
                auto dvdr = 48 * r2 * r6 * ( r6 - 0.5 ); 
                for (int d=0; d<3; ++d)
                    f[d] += dvdr * xr[d];
                // energy
                auto pot = 4 * r6 * ( r6 - 1 );
                en += pot;
            } // if

        } // for j
            
        //--- update acceleration and velocity    
        for (int d=0; d<3; ++d) {
            auto acc = f[d] * inv_mass;
            auto pos = 3*i + d;
            auto vd = v[pos] + 0.5 * dt * (acc + a[pos]);
            v[pos] = vd;
            a[pos] = acc;
            vsum[d] = vsum[d] + vd;
        }
            
    } // for i
    
    //--- remove linear momentum

    auto ninv = static_cast<real_t>(1) / size;
    for (int d=0; d<3; ++d)
        vsum[d] = vsum[d] * ninv;
       
    real_t ke0 = 0; 
    for (int i=0; i<size; ++i) {
        for (int d=0; d<3; ++d) {
            auto pos = 3*i + d;
            auto vd = v[pos] - vsum[d];
            v[pos] = vd;
            ke0 = ke0 + vd*vd; 
        }
    }
    
    //--- the temperature scaling
    auto num_dofs = 3*(size - 1);
    auto T0 = mass * ke0 / (boltzmann*num_dofs);
    auto fac = std::sqrt(temperature / T0 );
    
    for (int i=0; i<size; ++i) {
        for (int d=0; d<3; ++d) {
            auto pos = 3*i + d;
            v[pos] = v[pos] * fac;
        }
    }
    
}
}
