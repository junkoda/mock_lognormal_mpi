#include "grid.h"
#include <cmath>
#include <cassert>

Grid::Grid(const int nc_, const double boxsize_) :
  nc(nc_), ncz(nc_/2 + 1), boxsize(boxsize_)
{
  //size_t ncz= 2*(nc/2 + 1);

  ptrdiff_t local_nx_, local_x0_;
  ptrdiff_t total_size=
    fftw_mpi_local_size_3d(nc, nc, nc/2+1, MPI_COMM_WORLD,
			    &local_nx_, &local_x0_);

  fk= fftw_alloc_complex(total_size);
  fx= reinterpret_cast<double*>(fk);
  
  local_nx= local_nx_;
  local_x0= local_x0_;
  
  plan_forward=
    fftw_mpi_plan_dft_r2c_3d(nc, nc, nc, fx, fk,
			     MPI_COMM_WORLD, FFTW_ESTIMATE);
  
  plan_inverse=
    fftw_mpi_plan_dft_c2r_3d(nc, nc, nc, fk, fx,
			     MPI_COMM_WORLD, FFTW_ESTIMATE);
}

void Grid::fft_forward()
{
  assert(mode ==  grid_mode_x);
  fftw_mpi_execute_dft_r2c(plan_forward, fx, fk);
  mode = grid_mode_k;

  // f^(k) = \int f(x) e^-ikx d3x = sum f(x) dx^3 
  
  const double fac= pow(boxsize/nc, 3.0);
  const size_t ngrid= static_cast<size_t>(local_nx)*nc*2*(nc/2 + 1);
  for(size_t i=0; i<ngrid; ++i)
    fx[i] *= fac;
}

void Grid::fft_inverse()
{
  assert(mode ==  grid_mode_k);
  fftw_mpi_execute_dft_c2r(plan_inverse, fk, fx);
  mode = grid_mode_x;

  // f(x) = \int f(x) e^ikx d3k/(2pi)^3 = 1/V sum f(k) e^ikx
  const double fac= 1.0/(boxsize*boxsize*boxsize);
  const size_t ngrid= static_cast<size_t>(local_nx)*nc*2*(nc/2 + 1);

  for(size_t i=0; i<ngrid; ++i)
    fx[i] *= fac;
}

void Grid::clear()
{
  const size_t nx= local_nx;
  const size_t n= nc;
  const size_t ncz= 2*(nc/2 + 1);
  
  for(size_t ix=0; ix<nx; ++ix) {
    for(size_t iy=0; iy<n; ++iy) {
      for(size_t iz=0; iz<ncz; ++iz) {
	size_t index= (ix*n + iy)*ncz + iz;
	fx[index]= 0.0;
      }
    }
  }
}
