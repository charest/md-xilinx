#include "args.hpp"
#include "input.hpp"
#include "md.hpp"

#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using real_t = double;

using clock_timer_t = std::chrono::steady_clock;


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

    // create the model
    auto model = std::make_unique<md_t>(inputs);

    // total particles
    int num_cells = inputs.dims[0]*inputs.dims[1]*inputs.dims[2];
    int num_parts = model->size_lattice() * num_cells;

    // Place particles on a lattice
    std::vector<real_t> x(3*num_parts);
    std::vector<real_t> v(num_parts*3);
    std::vector<real_t> a(3*num_parts); 

    auto T0 = model->init_particles(x.data(), v.data(), a.data(), num_parts);

    std::cout << "md>> Initial temperature: " << T0 << std::endl;

    int output_counter = 0;
    
    if (inputs.output_freq) {
        std::cout << "md>> outputing: " << output_counter << std::endl;
    	output(output_counter++, x.data(), v.data(), num_parts, max_digits);
    }

    // start integration
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

        auto T0 = model->move_particles(
            x.data(),
            v.data(),
            a.data(),
            inputs.dt,
            num_parts );

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
