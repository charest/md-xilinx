#include "args.hpp"
#include "input.hpp"

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using real_t = double;

using clock_timer_t = std::chrono::steady_clock;

///////////////////////////////////////////////////////////////////////////////
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

///////////////////////////////////////////////////////////////////////////////
/// Return size of FCC lattice.  Note, only the repeating portion is needed.
int lattice_size_fcc()
{ return 4; }

///////////////////////////////////////////////////////////////////////////////
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

///////////////////////////////////////////////////////////////////////////////
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

///////////////////////////////////////////////////////////////////////////////
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


///////////////////////////////////////////////////////////////////////////////
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

void scale(
    real_t * v,
    real_t * fact,
    int num_part)
{
    for (int i=0; i<num_part; ++i)
        for (int d=0, pos=i*3; d<3; ++d, ++pos)
            v[pos] *= fact[d];

}

///////////////////////////////////////////////////////////////////////////////
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

///////////////////////////////////////////////////////////////////////////////
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

///////////////////////////////////////////////////////////////////////////////
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


///////////////////////////////////////////////////////////////////////////////
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

///////////////////////////////////////////////////////////////////////////////
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
    
///////////////////////////////////////////////////////////////////////////////
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

///////////////////////////////////////////////////////////////////////////////
/// Output solution in CSV format
void output(int n, const real_t * x, const real_t * v, int num_parts, int max_digits)
{
	std::stringstream ss;
	ss << "atoms";
	ss << std::setw(max_digits) << std::setfill('0') << n;
	ss << ".csv";

	std::ofstream file(ss.str());

	file << std::scientific;
	file.precision(12);
		
    file << "x, y, z, vx, vy, vz" << std::endl;

	for (int i=0; i<num_parts; ++i) {
		file << x[3*i] << ", " << x[3*i+1] << ", " << x[3*i+2] << ", ";
		file << v[3*i] << ", " << v[3*i+1] << ", " << v[3*i+2] << std::endl;
	}

}

///////////////////////////////////////////////////////////////////////////////
/// How many digits in an integer
int num_digits(int n) {
	if (n)
		return std::floor(std::log10(std::abs((real_t) n)) + 1);
	else
		return 1;
}

///////////////////////////////////////////////////////////////////////////////
/// Print usae
///////////////////////////////////////////////////////////////////////////////
void print_usage(char * cmd)
{
    std::cout << cmd << " [-h] [-f <xlcbin>] -i <input_file>" << std::endl;
    std::cout << "\t -h \t\t Print this help message." << std::endl;
    std::cout << "\t -f <xlcbin> \t Load <xlcbin> kernel image." << std::endl;
}


///////////////////////////////////////////////////////////////////////////////
/// Driver
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
    
    // command line arguments
    if(cmd_option_exists(argv, argv+argc, "-h"))
    {
        print_usage(argv[0]);
        return 0;
    }

    char * input_filename = get_cmd_option(argv, argv + argc, "-i");
    if (!input_filename) {
        std::cout << "No input specified." << std::endl;
        print_usage(argv[0]);
        return 1;
    }
   
#ifdef HAVE_VITIS 
    char * xlcbin_filename = get_cmd_option(argv, argv + argc, "-f");
#endif
    
    // Inputs
    std::cout << "md>> Loading inputs: '" << input_filename << "'" << std::endl;
    input_t inputs(input_filename);
    std::cout << std::endl << inputs << std::endl << std::endl;
    
    // get max digits expected
    int max_digits = num_digits(inputs.num_steps);

    // count total cells
    int num_cells = inputs.dims[0]*inputs.dims[1]*inputs.dims[2];

    // box size
    real_t box_length[] = {
        inputs.dims[0] * inputs.unit[0],
        inputs.dims[1] * inputs.unit[1],
        inputs.dims[2] * inputs.unit[2],
    };

    // Place particles on a lattice

    // total particles
    auto size_lattice = lattice_size_fcc();
    int num_parts = size_lattice * num_cells;

    // set positions
    std::vector<real_t> x(3*num_parts);
    init_fcc(x.data(), inputs.dims);

    // scale
    scale(x.data(), inputs.unit, num_parts);

    // Give particles random velocities
    std::vector<real_t> v(num_parts*3);
    init_random(v.data(), num_parts);

    // remove linear momentum
    linear_mom(v.data(), num_parts);

    // scale velocities
    auto T0 = scale_velocity(v.data(), inputs.mass, inputs.temperature,
                             inputs.boltzmann, num_parts); 

    std::cout << "md>> Initial temperature: " << T0 << std::endl;

    int output_counter = 0;
    
    if (inputs.output_freq) {
        std::cout << "md>> outputing: " << output_counter << std::endl;
    	output(output_counter++, x.data(), v.data(), num_parts, max_digits);
    }

    // start integration
    std::vector<real_t> a(3*num_parts, 0); 
    std::vector<real_t> f(3*num_parts);

    std::cout << "md>> " << std::string(73, '=') << std::endl;
    std::cout << "md>> " << "| " << std::setw(15) << "Iteration" << " | ";
    std::cout << std::setw(15) << "Soln Time" << " | ";
    std::cout << std::setw(15) << "Temperature" << " | ";  
    std::cout << std::setw(15) << "Elapsed Time, s" << " |";
    std::cout << std::endl;
    std::cout << "md>> " << std::string(73, '=') << std::endl;

#ifdef HAVE_VITIS
    // move everything over
#endif
	

    auto t = inputs.start_time;
    auto elapsed_start = clock_timer_t::now();

    int iter = 0;

    for (; iter<inputs.num_steps; ++iter, t+=inputs.dt)
    {

        // get new atom positions
        positions(v.data(), a.data(), x.data(), inputs.dt, num_parts);
        periodic(x.data(), box_length, num_parts);
        
        // compute forces
        auto en = forces(x.data(), f.data(), box_length, num_parts);

        // accleration
        accel(f.data(), v.data(), a.data(), inputs.mass, inputs.dt, num_parts);
        
        // remove linear momentum
        linear_mom(v.data(), num_parts);
        
        // scale velocities
        auto T0 = scale_velocity(
            v.data(), inputs.mass, inputs.temperature,
            inputs.boltzmann, num_parts); 

        if (inputs.stats_freq && ((iter+1) % inputs.stats_freq == 0)) {
            std::chrono::duration<real_t> elapsed = clock_timer_t::now() - elapsed_start;
            auto ss = std::cout.precision();
    	    std::cout << "md>> " << "| " << std::setw(15) << std::fixed << iter+1 << " | ";
    	    std::cout << std::setw(15) << std::scientific << t+inputs.dt << " | ";
    	    std::cout << std::setw(15) << std::scientific << T0 << " | ";  
    	    std::cout << std::setw(15) << std::scientific << elapsed.count() << " |";
    	    std::cout << std::endl; 
            std::cout.precision(ss);
            std::cout.unsetf( std::ios::scientific );
        }

        // output
        if (inputs.output_freq && ((iter+1) % inputs.output_freq == 0)) {
            std::cout << "md>> outputing: " << output_counter << std::endl;
            output(output_counter++, x.data(), v.data(), num_parts, max_digits);
        }
        
    } // nstep
	
    // output
    if (inputs.output_freq && (iter % inputs.output_freq != 0)) {
        std::cout << "md>> outputing: " << output_counter << std::endl;
        output(output_counter++, x.data(), v.data(), num_parts, max_digits);
    }

    std::chrono::duration<real_t> elapsed = clock_timer_t::now() - elapsed_start;
    std::cout << "md>> Total elapsed time (s): " << elapsed.count() << std::endl;

    return 0;
}
