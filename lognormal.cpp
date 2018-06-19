#include <iostream> // DEBUG
#include <vector>
#include <cmath>
#include <cassert>

#include <gsl/gsl_rng.h>

#include "comm.h" // DEBUG
#include "lognormal.h"

using namespace std;

namespace {
  void create_seedtable(const unsigned long seed, const int nc,
			vector<unsigned int>& seedtable);
} // unnamed namespace

void lognormal_compute_gaussian_power(Grid* const grid);


Grid* lognormal_create_grid(InputPower const * const ps,
			    const int nc, const double boxsize,
			    const double z, const double b)
{
  //
  // Create a grid of random lognormal field delta(x)
  //
  Grid* const grid= lognormal_create_power_grid(ps, nc, boxsize, 1.0);

  return grid;
}

Grid* lognormal_create_power_grid(InputPower const * const ps,
				  const int nc, const double boxsize,
				  const double pk_fac)
{
  //
  // Create a 3D grid of P_input(k)
  //
  Grid* const grid= new Grid(nc, boxsize);
  grid->clear();

  const int nx= grid->local_nx;
  const int nckz= nc/2 + 1;
  const double fac= 2.0*M_PI/boxsize;

  fftw_complex* const fk= grid->fk;
  
  for(int ix_local= 0; ix_local<nx; ++ix_local) {
    int ix= grid->local_x0 + ix_local;
    if(ix == nc/2) continue;
    int ikx= ix < nc/2 ? ix : ix - nc;
    for(int iy=0; iy<nc; ++iy) {
      if(iy == nc/2) continue;
      int iky= iy < nc/2 ? iy : iy - nc;
      size_t index_xy= (static_cast<size_t>(ix_local)*nc + iy)*nckz;

      int iz0= ix == 0 && iy == 0;
      for(int iz=iz0; iz<nc/2; ++iz) {
	int ikz= iz;

	double k= fac*sqrt(static_cast<double>(ikx*ikx + iky*iky + ikz*ikz));
	double P= pk_fac*ps->P(k);

	size_t index= index_xy + iz;

	fk[index][0]= P;
	fk[index][1]= 0.0;

	if(comm_this_node() == 1) printf("%e %e\n", k, P);
      }
    }
  }
  
  grid->mode= grid_mode_k;
  
  return grid;
}

void lognormal_compute_gaussian_power(Grid* const grid)
{
  //
  // Convert P_input(k) -> P_gaussian(k)
  //
  
  // P(k) => xi(r)
  grid->fft_inverse();

  const size_t nc= grid->nc;
  const size_t ncz= 2*(nc/2 + 1);
  const size_t nx= grid->local_nx;
  double* const xi= grid->fx;
  
  for(size_t ix=0; ix<nx; ++ix) {
    for(size_t iy=0; iy<nc; ++iy) {
      for(size_t iz=0; iz<nc; ++iz) {
	size_t index= (ix*nc + iy)*ncz + iz;
	assert(1.0 + xi[index] > 0.0);
	xi[index]= log(1.0 + xi[index]);
      }
    }
  }

  // xi_gaussian(r) => P_gaussian(r)
  grid->fft_forward();
}

void lognormal_create_gaussian_delta_k(const unsigned long seed,
				       Grid* const grid)
{
  // Input: grid as P(k)
  // Ouput: grid as delta(k)

  assert(grid->mode == grid_mode_k);

  vector<unsigned int> seedtable;
  create_seedtable(seed, grid->nc, seedtable);

  const int nc= grid->nc;
  const int nx= grid->local_nx;
  const int ix0= grid->local_x0;
  const size_t nckz= nc/2 + 1;
  const double boxsize= grid->boxsize;
  const double vol= boxsize*boxsize*boxsize;
  fftw_complex* const fk= (fftw_complex*) grid->fx;
  
  // P(k) = 1/V <delta(k) delta^*(k)
  gsl_rng* rng= gsl_rng_alloc(gsl_rng_ranlxd1);

  size_t negative= 0;
  double P_min= 0.0;

  for(int ix=0; ix<nc; ++ix) {
    int iix= nc - ix;
    if(iix == nc)
      iix= 0;
    
    if(!((ix0 <= ix  && ix  < ix0 + nx) ||
	 (ix0 <= iix && iix < ix0 + nx)))
      continue;

    for(int iy=0; iy<nc; ++iy) {
      int iiy = nc - iy;
      if(iiy == nc) iiy = 0;
      
      gsl_rng_set(rng, seedtable[ix * nc + iy]);
      
      for(int iz=0; iz<nc/2; ++iz) {
	double ampl= 1.0;
	do
	  ampl= gsl_rng_uniform(rng);
	while(ampl == 0.0);
	
	ampl= -log(ampl);
	  
	double phase= gsl_rng_uniform(rng)*2.0*M_PI;
	
	if(ix == nc/2 || iy == nc/2 || iz == nc/ 2)
	  continue;
	if(ix == 0 && iy == 0 && iz == 0)
	  continue;

	size_t index;
	if(ix0 <= ix && ix < ix0 + nx) {
	  assert(ix >= ix0);
	  index= ((ix - ix0)*nc + iy)*nckz + iz;
	}
	else if(ix0 <= iix && iix < ix0 + nx) {
	  assert(iix - ix0);
	  index= ((iix - ix0)*nc + iy)*nckz + iz;
	}
	else {
	  assert(false);
	}


	double delta2= vol*fk[index][0];

	double delta_k_mag= 0.0;
	if(fk[index][0] < P_min)
	  P_min= fk[index][0];
	
	if(fk[index][0] > 0.0)
	  delta_k_mag= sqrt(ampl*delta2);
	else
	  negative++;
	
	if(iz > 0) {
	  if(ix0 <= ix && ix < ix0 + nx) {
	    assert(ix >= ix0);
	    fk[((ix - ix0)*nc + iy)*nckz + iz][0]= delta_k_mag*cos(phase);
	    fk[((ix - ix0)*nc + iy)*nckz + iz][1]= delta_k_mag*sin(phase);
	  }
	}
	else {
	  // iz=0 plane: assign also delta(-k) = delta(k)^* [reality condition]
	  if(ix == 0) {
	    if(iy >= nc/2) {
	      continue;
	    }
	    else {
	      if(ix0 <= ix && ix < ix0 + nx) {
		assert(ix >= ix0);
		
		fk[((ix - ix0)*nc + iy)*nckz + iz][0]= delta_k_mag*cos(phase);
		fk[((ix - ix0)*nc + iy)*nckz + iz][1]= delta_k_mag*sin(phase);
		
		fk[((ix - ix0)*nc + iiy)*nckz + iz][0]= delta_k_mag*cos(phase);
		fk[((ix - ix0)*nc + iiy)*nckz + iz][1]= -delta_k_mag*sin(phase);	      }
	    }
	  }
	  else {
	    if(ix >= nc/2)
	      continue;
	    else {

	      if(ix0 <= ix && ix < ix0 + nx) {
		fk[((ix - ix0)*nc + iy)*nckz + iz][0]= delta_k_mag*cos(phase);
		fk[((ix - ix0)*nc + iy)*nckz + iz][1]= delta_k_mag*sin(phase);
	      }

	      if(ix0 <= iix && iix < ix0 + nx) {
		fk[((iix - ix0)*nc + iiy)*nckz + iz][0]= delta_k_mag*cos(phase);
		fk[((iix - ix0)*nc + iiy)*nckz + iz][1]= -delta_k_mag*sin(phase);
	      }
	    }
	  }
	}
      }
    }
  }

  gsl_rng_free(rng);

  //fprintf(stderr, "P_min= %e\n", P_min);
  //fprintf(stderr, "negative P(k): %zu\n", negative);  
}

namespace {
  
void create_seedtable(const unsigned long seed, const int nc,
		      vector<unsigned int>& seedtable)
{
  seedtable.resize(nc*nc);
  
  //
  // Setup random seeds (From NGen-IC)
  //
  gsl_rng* random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(random_generator, seed);

  for(int i=0; i<nc/2; i++) {
    for(int j=0; j<i; j++)
      seedtable[i * nc + j] = 0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i+1; j++)
      seedtable[j * nc + i] = 0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i; j++)
      seedtable[(nc - 1 - i) * nc + j] = 
	0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i+1; j++)
      seedtable[(nc - 1 - j) * nc + i] = 
	0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i; j++)
      seedtable[i * nc + (nc - 1 - j)] = 
	0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i+1; j++)
      seedtable[j * nc + (nc - 1 - i)] = 
	0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i; j++)
      seedtable[(nc - 1 - i) * nc + (nc - 1 - j)] = 
	0x7fffffff * gsl_rng_uniform(random_generator);

    for(int j=0; j<i+1; j++)
      seedtable[(nc - 1 - j) * nc + (nc - 1 - i)] = 
	0x7fffffff * gsl_rng_uniform(random_generator);
  }
  
  gsl_rng_free(random_generator);
}

} // End of unnamed namespace
