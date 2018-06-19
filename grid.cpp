#include "grid.h"

Grid::Grid(const int nc_, const double boxsize_) :
  nc(nc_), boxsize(boxsize_)
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
  fftw_mpi_execute_dft_r2c(plan_forward, fx, fk);
}

void Grid::fft_inverse()
{
  fftw_mpi_execute_dft_c2r(plan_inverse, fk, fx);
}
