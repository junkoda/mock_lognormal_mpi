#include <iostream>
#include <cstdio>
#include <cassert>

#include "comm.h"
#include "input_power.h"
#include "grid.h"
#include "lognormal.h"
#include "hdf5_write.h"

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
    ("nc", value<int>()->default_value(256), "number of grids per dimension")
    ("boxsize", value<double>()->default_value(1000.0),
                "length of the box on a side")
    ("write-pk-input", value<string>()->default_value("p3d.h5"),
     "write 3D grid of input power spectrum")
    ("filename", value<string>(), "P(k) file name")
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

  InputPower* const ps= new InputPower(filename.c_str());

  Grid* const grid= lognormal_create_power_grid(ps, nc, boxsize, 1.0);

  
  if(vm.count("write-pk-input")) {
    hdf5_write_grid_complex(vm["write-pk-input"].as<string>().c_str(),
			    grid);
  }

  comm_finalise();
  
  return 0;
}
