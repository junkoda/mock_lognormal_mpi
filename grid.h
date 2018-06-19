#ifndef GRID_H
#define GRID_H 1

#include <fftw3-mpi.h>

enum FFTMode {grid_mode_unknown, grid_mode_x, grid_mode_k};

class Grid {
 public:
  Grid(const int nc_, const double boxsize);
  ~Grid();

  void fft_forward();
  void fft_inverse();
  void clear();

  double* fx;
  fftw_complex* fk;
  const int nc;
  const double boxsize;
  int local_nx, local_x0;
  
  FFTMode mode;
  double x0_box[3];

  
 private:
  fftw_plan plan_forward;
  fftw_plan plan_inverse;

};

#endif
