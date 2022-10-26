#include <iostream>
#include <random>
#include <vector>

using real_t = double;

class rand_t {
    std::mt19937 gen_;
    std::uniform_real_distribution<> dist_;


public:
    template<typename T>
    rand_t(const T & seed) :
      gen_(seed), dist_(-0.5, 0.5)
    {}

    real_t operator()()
    { return dist_(gen_); }
};

int lattice_size_fcc()
{ return 4; }

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
                for (int l=0, lpos=0; l<lattice_size; ++l, ++pos)
                {
                    x[3*pos + 0] = i + lattice[3*lpos+0];
                    x[3*pos + 1] = j + lattice[3*lpos+1];
                    x[3*pos + 2] = k + lattice[3*lpos+2];
                }

}

void init_random(
    real_t * v,
    int num_part)
{
    rand_t rand(0);

    for (int i=0; i<num_part; ++i)
    {
        for (int d=0; d<3; ++d) {
            auto pos = i*3 + d;
            v[pos] = rand();
        }
    }
}

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


void scale(
    real_t * v,
    real_t fact,
    int num_part)
{
    for (int i=0; i<num_part; ++i)
        for (int d=0, pos=i*3; d<3; ++d, ++pos)
            v[pos] *= fact;

}

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


real_t forces(
    const real_t * x,
    real_t * f,
    const real_t * box_length,
    int num_part)
{
    real_t en = 0;
    
    for (int pos=0; pos<num_part*3; ++pos)
        f[pos] = 0;

    for (int i=0; i<num_part; ++i)
        for (int j=0; j<num_part; ++j)
        {
            if (i != j) {
                // distance
                real_t xr[3] = {
                    x[3*i + 0] - x[3*j + 0],
                    x[3*i + 1] - x[3*j + 1],
                    x[3*i + 2] - x[3*j + 2],
                };
                // periodic boundary conditions
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
                    f[i*3 + d] += dvdr * xr[d];
                // energy
                auto pot = 4 * r6 * ( r6 - 1 );
                en += pot;
            }
        }

    return en / 2;
}

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

int main(int argc, char* argv[])
{
    // 0. Inputs

    // temperature
    real_t temperature = 10;
    // boltzmann constant
    real_t boltzmann = 1;
    // mass
    real_t mass = 1;
    // grid
    int dims[] = {10, 10, 10};
    // lattice constant
    real_t lattice_constant = std::pow( 2., 2./3. );
    // time step
    real_t dt = 0.001;
    // number of steps
    int num_steps = 1;

    // count total cells
    int num_cells = dims[0]*dims[1]*dims[2];

    // box size
    real_t box_length[] = {
        dims[0] * lattice_constant,
        dims[1] * lattice_constant,
        dims[2] * lattice_constant,
    };
    
    // 1. Place particles on a lattice

    // total particles
    auto size_lattice = lattice_size_fcc();
    int num_parts = size_lattice * num_cells;

    // set positions
    std::vector<real_t> x(3*num_parts);
    init_fcc(x.data(), dims);

    // scale
    for (auto & xi : x )
        xi *= lattice_constant;

    // 2. Give particles random velocities
    std::vector<real_t> v(num_parts);
    init_random(v.data(), num_parts);

    // remove linear momentum
    linear_mom(v.data(), num_parts);

    // scale velocities
    auto num_dofs = 3*(num_parts - 1);
    auto ke0 = dot_product(v.data(), num_parts);
    auto T0 = mass * ke0 / (boltzmann*num_dofs);
    auto fac = std::sqrt(temperature / T0 );
    scale(v.data(), fac, num_parts);

    // start integration
    std::vector<real_t> a(3*num_parts, 0); 
    for (int n=0; n<num_steps; ++n)
    {

        // 3. get new atom positions
        positions(v.data(), a.data(), x.data(), dt, num_parts);
        periodic(x.data(), box_length, num_parts);

        // 4. compute forces
        std::vector<real_t> f(3, num_parts);
        auto en = forces(x.data(), f.data(), box_length, num_parts);

        // 5. accleration
        accel(f.data(), v.data(), a.data(), mass, dt, num_parts);

        // remove linear momentum
        linear_mom(v.data(), num_parts);

        // scale velocities
        auto num_dofs = 3*(num_parts - 1);
        auto ke0 = dot_product(v.data(), num_parts);
        auto T0 = mass * ke0 / (boltzmann*num_dofs);
        auto fac = std::sqrt(temperature / T0 );
        scale(v.data(), fac, num_parts);

        std::cout << "Temperature = " << T0 << std::endl;

    } // nstep

    return 0;
}