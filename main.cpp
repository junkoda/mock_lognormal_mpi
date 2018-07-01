#include <iostream>
#include <cstdio>
#include <cassert>
#include <gsl/gsl_rng.h>

#include "comm.h"
#include "msg.h"
#include "input_power.h"
#include "grid.h"
#include "lognormal.h"
#include "config.h"
#include "hdf5_write.h"
#include "mass_assignment.h"

#include <boost/program_options.hpp>

using namespace std;
using namespace boost::program_options;


int main(int argc, char* argv[])
{
  comm_init(&argc, &argv);

  options_description
    opt("mpirun -n 1 mock_lognormal_mpi [options] <P(k) filename>");
  
  opt.add_options()
    ("help,h", "display this help")
    ("filename", value<string>(), "P(k) file name")
    ("nc", value<int>()->default_value(256), "number of grids per dimension")
    ("boxsize", value<double>()->default_value(1000.0),
                "length of the box on a side")
    ("seed", value<unsigned long>()->default_value(1), "random seed")
    ("nbar-uniform", value<double>()->default_value(1.0e-4),
     "number density of the uniform mock")
    ("mas", value<int>()->default_value(2), "mass assignment scheme; 1=NGP, 2=CIC, 3=TSC")
    ("write-pk-input", value<string>(),
       "write 3D grid of input power spectrum")
    ("write-delta-k", value<string>(), "write 3D grid of delta(k)")
    ("write-lognormal-density", value<string>(), "write 3D grid of lognormal n(x)")
    ("write-mock-particles", value<string>(),
     "write mock particle position and velocity")
    ("write-mock-density", value<string>(), "write 3D grid of mock n(x)")
    ;

  positional_options_description p;
  p.add("filename", -1);

  variables_map vm;
  store(command_line_parser(argc, argv).options(opt).positional(p).run(), vm);
  notify(vm);

  if(vm.count("help") || ! vm.count("filename")) {
    cout << opt;
    return 0;
  }

  const string filename = vm["filename"].as<string>();
  const int nc= vm["nc"].as<int>(); assert(nc > 0);
  const double boxsize= vm["boxsize"].as<double>(); assert(boxsize > 0.0);

  //
  // Create 3D grid P(k)
  //
  InputPower* const ps= new InputPower(filename.c_str());

  Grid* const grid= lognormal_create_power_grid(ps, nc, boxsize, 1.0);

  
  if(vm.count("write-pk-input")) {
    hdf5_write_grid_complex(vm["write-pk-input"].as<string>().c_str(),
			    grid);
  }

  //
  // Convert to P_guassian(k)
  //
  //lognormal_compute_gaussian_power(grid);
  lognormal_convert_to_gaussian_power(grid);
  
  //
  // Convert P_gaussian(k) to random realisation delta_k(k)
  //
  const unsigned long seed= vm["seed"].as<unsigned long>();
  lognormal_generate_gaussian_delta_k(seed, grid);

  
  if(vm.count("write-delta-k")) {
    hdf5_write_grid_complex(vm["write-delta-k"].as<string>().c_str(),
			    grid);
  }

  //
  // Convert Gaussian density to lognormal density
  // 1 + delta_gaussian(x) => 1 + delta_lognormal(x)
  //
  msg_printf(msg_debug, "FFT to gaussian density n_gaussian(x)\n");

  grid->fft_inverse();
  
  lognormal_convert_to_lognormal_density(grid);

  if(vm.count("write-lognormal-density")) {
    hdf5_write_grid_real(vm["write-lognormal-density"].as<string>().c_str(),
			 grid);
  }


  //
  // TODO: create velocity grids
  //

  //
  // Generate uniform mock
  //
  msg_printf(msg_debug, "generating mock particles\n");
  
  //gsl_rng* rng= gsl_rng_alloc(gsl_rng_ranlxd1);
  //gsl_rng_set(rng, 40001 + seed + 105*comm_this_node());
    
  const double nbar= vm["nbar-uniform"].as<double>();

  Particles particles;
  lognormal_generate_particles_periodic(nbar, seed, grid,
					0, 0, 0, particles);

  cerr << "comm_sum\n";

  long long np_mock= particles.size();
  np_mock= comm_sum(np_mock);
  msg_printf(msg_info, "%lld particles generated in the periodic box\n",
	     np_mock);

  if(vm.count("write-mock-particles")) {
    hdf5_write_grid_particles(vm["write-mock-particles"].as<string>().c_str(),
			      particles);
  }


  /*
  //
  // Compute box power spectrum
  //
  Grid* const grid_ps= new Grid(nc/2, boxsize);
  grid_ps->clear();
  
  msg_printf(msg_debug, "computing real-space power spectrum\n");

  const int mas= vm["mas"].as<int>();
  mass_assignment(particles, mas, grid_ps, 0);

  if(vm.count("write-mock-density")) {
    hdf5_write_grid_real(vm["write-mock-density"].as<string>().c_str(),
			 grid_ps);
  }
  */  


  
  

  comm_finalise();
  
  return 0;
}
