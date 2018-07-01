#include <cmath>

#include "msg.h"
#include "comm.h"
#include "mass_assignment.h"


namespace {
//
// Mass assignment functions
//  x: position relative to the cubic box corner in units of grid spacing;
//     i.e., (0, 0, 0) and (nc, nc, nc) is the edge of the FFT grid,
//     where nc is the number of grids per dimension

struct NGP {
  void operator()(const double x[], Grid* const grid) const {
    int ix[3];
    for(int i=0; i<3; ++i) {
      ix[i] = (int) floor(x[i] + 0.5);
      ix[i] = (ix[i] + grid->nc) % grid->nc;
    }
    

    if(grid->local_x0 <= ix[0] && ix[0] < grid->local_x0 + grid->local_nx)
      grid->add(ix[0], ix[1], ix[2], 1.0);
  }
  static const int n_mas = 1;
};
    

struct CIC {
  void operator()(const double x[], Grid* const grid) const {
    int ix[3], ix0[3], ix1[3];
    double w0[3], w1[3];

    for(int k=0; k<3; ++k) {
      ix[k] = (int) floor(x[k]);
      ix0[k]= (ix[k] + grid->nc) % grid->nc;    // left grid point (periodic)
      ix1[k]= (ix[k] + 1 + grid->nc) % grid->nc;// right grid point (periodic)

      w1[k] = x[k] - ix[k];              // CIC weight to right grid point
      w0[k] = 1.0 - w1[k];                 //               left grid point
    }

    if(grid->local_x0 <= ix0[0] && ix0[0] < grid->local_x0 + grid->local_nx) {
      grid->add(ix0[0], ix0[1], ix0[2], w0[0]*w0[1]*w0[2]);
      grid->add(ix0[0], ix1[1], ix0[2], w0[0]*w1[1]*w0[2]);
      grid->add(ix0[0], ix0[1], ix1[2], w0[0]*w0[1]*w1[2]);
      grid->add(ix0[0], ix1[1], ix1[2], w0[0]*w1[1]*w1[2]);
    }

    if(grid->local_x0 <= ix1[0] && ix1[0] < grid->local_x0 + grid->local_nx) {
      grid->add(ix1[0], ix0[1], ix0[2], w1[0]*w0[1]*w0[2]);
      grid->add(ix1[0], ix1[1], ix0[2], w1[0]*w1[1]*w0[2]);
      grid->add(ix1[0], ix0[1], ix1[2], w1[0]*w0[1]*w1[2]);
      grid->add(ix1[0], ix1[1], ix1[2], w1[0]*w1[1]*w1[2]);
    }
  }

  static const int n_mas = 2;
};

struct TSC {
  void operator()(const double x[], Grid* const grid) const {
    const int ix_left= grid->local_x0;
    const int ix_right= grid->local_x0 + grid->local_nx;

    int ix[3], ix0[3], ix1[3], ix2[3];
    double w0[3], w1[3], w2[3];

    for(int k=0; k<3; ++k) {
      ix[k] = (int) floor(x[k] + 0.5);
      ix0[k]= (ix[k] - 1 + grid->nc) % grid->nc;
      ix1[k]= (ix[k]     + grid->nc) % grid->nc;
      ix2[k]= (ix[k] + 1 + grid->nc) % grid->nc;

      double dx1 = x[k] - ix[k];
      double dx2 = 0.5 - dx1;
      
      w0[k] = 0.5*dx2*dx2;
      w1[k] = 0.75 - dx1*dx1;
      w2[k] = 0.25 - 0.5*dx2*dx2 + dx1*dx1;
    }

    if(ix_left <= ix0[0] && ix0[0] < ix_right) {
      grid->add(ix0[0], ix0[1], ix0[2], w0[0]*w0[1]*w0[2]);
      grid->add(ix0[0], ix0[1], ix1[2], w0[0]*w0[1]*w1[2]);
      grid->add(ix0[0], ix0[1], ix2[2], w0[0]*w0[1]*w2[2]);
      grid->add(ix0[0], ix1[1], ix0[2], w0[0]*w1[1]*w0[2]);
      grid->add(ix0[0], ix1[1], ix1[2], w0[0]*w1[1]*w1[2]);
      grid->add(ix0[0], ix1[1], ix2[2], w0[0]*w1[1]*w2[2]);
      grid->add(ix0[0], ix2[1], ix0[2], w0[0]*w2[1]*w0[2]);
      grid->add(ix0[0], ix2[1], ix1[2], w0[0]*w2[1]*w1[2]);
      grid->add(ix0[0], ix2[1], ix2[2], w0[0]*w2[1]*w2[2]);
    }
    if(ix_left <= ix1[0] && ix1[0] < ix_right) {
      grid->add(ix1[0], ix0[1], ix0[2], w1[0]*w0[1]*w0[2]);
      grid->add(ix1[0], ix0[1], ix1[2], w1[0]*w0[1]*w1[2]);
      grid->add(ix1[0], ix0[1], ix2[2], w1[0]*w0[1]*w2[2]);
      grid->add(ix1[0], ix1[1], ix0[2], w1[0]*w1[1]*w0[2]);
      grid->add(ix1[0], ix1[1], ix1[2], w1[0]*w1[1]*w1[2]);
      grid->add(ix1[0], ix1[1], ix2[2], w1[0]*w1[1]*w2[2]);
      grid->add(ix1[0], ix2[1], ix0[2], w1[0]*w2[1]*w0[2]);
      grid->add(ix1[0], ix2[1], ix1[2], w1[0]*w2[1]*w1[2]);
      grid->add(ix1[0], ix2[1], ix2[2], w1[0]*w2[1]*w2[2]);
    }
    if(ix_left <= ix2[0] && ix2[0] < ix_right) {
      grid->add(ix2[0], ix0[1], ix0[2], w2[0]*w0[1]*w0[2]);
      grid->add(ix2[0], ix0[1], ix1[2], w2[0]*w0[1]*w1[2]);
      grid->add(ix2[0], ix0[1], ix2[2], w2[0]*w0[1]*w2[2]);
      grid->add(ix2[0], ix1[1], ix0[2], w2[0]*w1[1]*w0[2]);
      grid->add(ix2[0], ix1[1], ix1[2], w2[0]*w1[1]*w1[2]);
      grid->add(ix2[0], ix1[1], ix2[2], w2[0]*w1[1]*w2[2]);
      grid->add(ix2[0], ix2[1], ix0[2], w2[0]*w2[1]*w0[2]);
      grid->add(ix2[0], ix2[1], ix1[2], w2[0]*w2[1]*w1[2]);
      grid->add(ix2[0], ix2[1], ix2[2], w2[0]*w2[1]*w2[2]);
    }
  }

  static const int n_mas = 3;
};

template<typename MAS>
void assign_template(const Particles& v,
		     MAS f,
		     const double offset,
		     Grid* const grid)
{
  assert(grid->boxsize > 0);
  assert(grid->nc > 0);

  const int nc= grid->nc;
  const double boxsize= grid->boxsize;
  const double dx_inv= nc/boxsize;
  //const double dx= boxsize/nc;
  //const double x0= offset;

  for(size_t i=0; i<v.size(); ++i) {
    double rx[3];
    
    rx[0] = v[i].x[0]*dx_inv - offset;
    rx[1] = v[i].x[1]*dx_inv - offset;
    rx[2] = v[i].x[2]*dx_inv - offset;
    
    f(rx, grid);
  }

  long long np= v.size();
  np= comm_sum(np);

  grid->np += np;
  grid->n_mas = f.n_mas;
}

void assign(const Particles& v, const int mas, const double offset,
	    Grid* const grid)
{
  switch(mas) {
  case 1:
    assign_template(v, NGP(), offset, grid);
    break;

  case 2:
    assign_template(v, CIC(), offset, grid);
    break;

  case 3:
    assign_template(v, TSC(), offset, grid);
    break;

  default:
    msg_abort("Error: unknown mass assignment scheme %d\n", mas);

  }
}
  
} // unnamed namespace

void mass_assignment(const Particles& v,
		     const int mas,
		     Grid* const grid,
		     Grid* const grid_shifted)
{
  // 1. Mass assignment for local particles
  assign(v, mas, 0.5, grid);

  // 2. Mass assignment for particles in other MPI nodes
  // TODO




  grid->mode= grid_mode_x;
  if(grid_shifted)
    grid_shifted->mode= grid_mode_x;
}
