#ifndef GRID_H
#define GRID_H 1

#include <cassert>
#include <fftw3-mpi.h>


enum FFTMode {grid_mode_unknown, grid_mode_x, grid_mode_k};

class Grid {
 public:
  Grid(const int nc_, const double boxsize);
  ~Grid();

  void fft_forward();
  void fft_inverse();
  void clear();

  void add(const int ix, const int iy, const int iz,
	   const double w) {
    int ix_local= ix - local_x0;
    assert(0 <= ix_local && ix_local < local_nx); //DEBUG
    size_t index= (ix*static_cast<size_t>(local_nx) + iy)*ncz;
    fx[index] += w;
  }

  double* fx;
  fftw_complex* fk;
  const int nc, ncz;
  const double boxsize;
  int local_nx, local_x0;
  
  FFTMode mode;
  //double x0_box[3];

  size_t np;
  int n_mas;
  
 private:
  fftw_plan plan_forward;
  fftw_plan plan_inverse;

};

#endif
