#include "md.hpp"

#include "input.hpp"
#include "random.hpp"

#include <cmath>
#include <vector>

/// Return size of FCC lattice.  Note, only the repeating portion is needed.
int lattice_size_fcc()
{ return 4; }

/// Place particles on the FCC lattice
void init_fcc(
    real_t * x,
    const int * num_cells)
{
    auto lattice_size = lattice_size_fcc();

    std::vector<real_t> lattice;
    lattice.reserve(3*lattice_size);

    lattice.push_back( 0 );
    lattice.push_back( 0 );
    lattice.push_back( 0 );

    lattice.push_back( 0.5 );
    lattice.push_back( 0.5 );
    lattice.push_back( 0 );

    lattice.push_back( 0 );
    lattice.push_back( 0.5 );
    lattice.push_back( 0.5 );

    lattice.push_back( 0.5 );
    lattice.push_back( 0 );
    lattice.push_back( 0.5 );

    for (int i=0, pos=0; i<num_cells[0]; ++i)
        for (int j=0; j<num_cells[1]; ++j)
            for (int k=0; k<num_cells[2]; ++k)
                for (int l=0; l<lattice_size; ++l, ++pos)
                {
                    x[3*pos + 0] = i + lattice[3*l+0];
                    x[3*pos + 1] = j + lattice[3*l+1];
                    x[3*pos + 2] = k + lattice[3*l+2];
                }

}

/// Initialize an array with random numbers between -0.5 and 0.5
void init_random(
    real_t * v,
    int num_part)
{
    rand_t rand(-0.5, 0.5, 0);

    for (int i=0, pos=0; i<num_part; ++i)
        for (int d=0; d<3; ++d, ++pos)
            v[pos] = rand();
}

/// Compute a dot product of an array with itself
real_t dot_product(
    real_t * v,
    int num_part)
{
    real_t sumv2 = 0;

    for (int i=0; i<num_part; ++i)
    {
        for (int d=0, pos=i*3; d<3; ++d, ++pos) {
            sumv2 += v[pos]*v[pos]; // kinetic energy
        }
    }

    return sumv2;

}

/// Scale an array
void scale(
    real_t * v,
    real_t fact,
    int num_part)
{
    for (int i=0; i<num_part; ++i)
        for (int d=0, pos=i*3; d<3; ++d, ++pos)
            v[pos] *= fact;

}

/// Scale an array
void scale(
    real_t * v,
    const real_t * fact,
    int num_part)
{
    for (int i=0; i<num_part; ++i)
        for (int d=0, pos=i*3; d<3; ++d, ++pos)
            v[pos] *= fact[d];

}

/// Move particles
void positions(
    const real_t * v,
    const real_t * a,
    real_t * x,
    real_t dt,
    int num_part)
{  
    auto dt_sqr = dt*dt;
    for (int pos=0; pos<num_part*3; ++pos)
        x[pos] += v[pos]*dt + 0.5*dt_sqr*a[pos];
}

/// Apply periodic boundaries
void periodic(real_t * x, const real_t * box_length, int num_part)
{
    for (int p=0, pos=0; p<num_part; ++p)
    {
        for (int d=0; d<3; ++d, ++pos)
        {
            auto & xi = x[pos];
            auto len = box_length[d];
            if (xi > len)
                xi -= len;
            else if (xi < 0)
                xi += len;
        }
    }

}


/// Compute forces
real_t forces(
    const real_t * x,
    real_t * f,
    const real_t * box_length,
    int num_part)
{
    real_t en = 0;
    
    for (int pos=0; pos<num_part*3; ++pos)
        f[pos] = 0;

    for (int i=0; i<num_part-1; ++i)
        for (int j=i+1; j<num_part; ++j)
        {
            //if (i != j) {
                // distance
                real_t xr[3] = {
                    x[3*i + 0] - x[3*j + 0],
                    x[3*i + 1] - x[3*j + 1],
                    x[3*i + 2] - x[3*j + 2],
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
                for (int d=0; d<3; ++d) {
                    f[i*3 + d] += dvdr * xr[d];
                    f[j*3 + d] -= dvdr * xr[d];
                }
                // energy
                auto pot = 4 * r6 * ( r6 - 1 );
                en += pot;
            //}
        }

    return en / 2;
}

/// Update acceleration and velocity
void accel(
    const real_t * f,
    real_t * v,
    real_t * a,
    real_t mass,
    real_t dt,
    int num_parts)
{
    auto inv_mass = 1 / mass;
    for (int pos=0; pos<3*num_parts; ++pos) {
        auto acc = f[pos] * inv_mass;
        v[pos] += 0.5 * dt * (acc + a[pos]);
        a[pos] = acc;
    }
}

/// Remove linear momentum
void linear_mom(real_t * v, int num_parts)
{

    real_t vsum[] = {0, 0, 0};
    for (int p=0, pos=0; p<num_parts; ++p)
        for (int d=0; d<3; ++d, ++pos)
            vsum[d] += v[pos];
    
    for (int d=0; d<3; ++d)
        vsum[d] /= num_parts;
    
    for (int p=0, pos=0; p<num_parts; ++p)
        for (int d=0; d<3; ++d, ++pos)
            v[pos] -= vsum[d];
    
}
    
/// Scale velocityies
real_t scale_velocity(
    real_t * v,
    real_t mass,
    real_t temperature,
    real_t boltzmann,
    int num_parts)
{
    auto num_dofs = 3*(num_parts - 1);
    auto ke0 = dot_product(v, num_parts);
    auto T0 = mass * ke0 / (boltzmann*num_dofs);
    auto fac = std::sqrt(temperature / T0 );
    scale(v, fac, num_parts);
    return T0;
}

md_t::md_t(const input_t & in)
{
    for (int d=0; d<3; ++d) {
        dims_[d] = in.dims[d];
        unit_[d] = in.unit[d];
        box_length_[d] = in.dims[d] * in.unit[d];
    }
    temperature_ = in.temperature;
    boltzmann_ = in.boltzmann;
    mass_ = in.mass;
}

int md_t::size_lattice() const
{ return lattice_size_fcc(); }

real_t md_t::init_particles(
    real_t * x,
    real_t * v,
    real_t * a,
    int num_parts) const
{
    // set positions
    init_fcc(x, dims_);

    // scale
    scale(x, unit_, num_parts);

    // Give particles random velocities
    init_random(v, num_parts);

    // remove linear momentum
    linear_mom(v, num_parts);

    // scale velocities
    auto T0 = scale_velocity(v, mass_, temperature_, boltzmann_, num_parts); 

    std::fill_n(a, 3*num_parts, 0);

    return T0;

}

real_t md_t::move_particles(
    real_t * x,
    real_t * v,
    real_t * a,
    real_t dt,
    int num_parts) const
{
    std::vector<real_t> f(3*num_parts);

    // get new atom positions
    positions(v, a, x, dt, num_parts);
    periodic(x, box_length_, num_parts);

    // compute forces
    auto en = forces(x, f.data(), box_length_, num_parts);

    // accleration
    accel(f.data(), v, a, mass_, dt, num_parts);

    // remove linear momentum
    linear_mom(v, num_parts);

    // scale velocities
    auto T0 = scale_velocity(
        v, mass_, temperature_,
        boltzmann_, num_parts);
    return T0;
}
