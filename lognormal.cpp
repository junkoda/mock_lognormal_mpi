#include <iostream> // DEBUG
#include <vector>
#include <cmath>
#include <cassert>

//#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "comm.h" // DEBUG
#include "msg.h"
#include "growth.h"
#include "lognormal.h"

using namespace std;

namespace {
  void create_seedtable(const unsigned long seed, const int nc,
			vector<unsigned int>& seedtable);

} // unnamed namespace

// inline function
static inline void get_random_phase(gsl_rng* rng,
				    double* ampl_out, double* phase_out) {
  // Generate random amplitude and phase of random Gaussian field
  double ampl;
  do
    ampl= gsl_rng_uniform(rng);
  while(ampl == 0.0);
	
  *ampl_out= -log(ampl);
  *phase_out= gsl_rng_uniform(rng)*2.0*M_PI;
}


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

	assert(P > 0.0);
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
  // using xi_gaussian(r)= log(1.0 + xi_input(r))
  
  // P(k) => xi(r)
  grid->fft_inverse();

  const size_t nc= grid->nc;
  const size_t ncz= 2*(nc/2 + 1);
  const size_t nx= grid->local_nx;
  double* const xi= grid->fx;
  
  for(size_t ix_local=0; ix_local<nx; ++ix_local) {
    for(size_t iy=0; iy<nc; ++iy) {
      for(size_t iz=0; iz<nc; ++iz) {
	size_t index= (ix_local*nc + iy)*ncz + iz;
	assert(1.0 + xi[index] > 0.0);
	xi[index]= log(1.0 + xi[index]);
      }
    }
  }

  // xi_gaussian(r) => P_gaussian(k)
  grid->fft_forward();
}


void lognormal_create_gaussian_delta_k(const unsigned long seed,
				       Grid* const grid)
{
  // Convert the grid of P(k) to a random Gaussian realisation delta(k)
  
  // Input: grid as P(k)
  //        seed: random seed
  //
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

  //
  // iz > 0
  //   all modes are independent
  for(int ix_local=0; ix_local<nx; ++ix_local) {
    int ix= ix0 + ix_local;
    if(ix == nc/2) continue;
    for(int iy=0; iy<nc; ++iy) {
      if(iy == nc/2) continue;
      
      gsl_rng_set(rng, seedtable[ix*nc + iy]);

      double ampl, phase;
      get_random_phase(rng, &ampl, &phase); // random phase for iz=0
      
      for(int iz=1; iz<nc/2; ++iz) {
	get_random_phase(rng, &ampl, &phase);

	size_t index= (ix_local*nc + iy)*nckz + iz;
	double delta2= vol*fk[index][0];

	double delta_k_mag= 0.0;
	if(delta2 > 0.0)
	  delta_k_mag= sqrt(ampl*delta2);
	else
	  negative++;

	fk[index][0]= delta_k_mag*cos(phase);
	fk[index][1]= delta_k_mag*sin(phase);
      }
    }
  }

  //
  // iz = 0
  //   Half of the modes are independent
  //   reality condition delta(-k) = delta(k)^* needs to be satisfied

  for(int ix=0; ix<nc/2; ++ix) {
    int ix_local= ix - ix0;
    int iix= nc - ix;
    int iix_local= iix - ix0;
    
    for(int iy=0; iy<nc; ++iy) {
      if(iy == nc/2) continue;
      else if(ix == 0 && iy > nc/2) continue;
      else if(ix == 0 && iy == 0) continue;

      // all modes here are independent

      gsl_rng_set(rng, seedtable[ix*nc + iy]);

      double ampl, phase;
      get_random_phase(rng, &ampl, &phase);

      if(0 <= ix_local && ix_local < nx) {
	size_t index= (ix_local*nc + iy)*nckz; // iz=0
	double delta2= vol*fk[index][0];
	double delta_k_mag= 0.0;
	if(delta2 > 0.0)
	  delta_k_mag= sqrt(ampl*delta2);
	else
	  negative++;

	fk[index][0]= delta_k_mag*cos(phase);
	fk[index][1]= delta_k_mag*sin(phase);
      }
      
      if(0 <= iix_local && iix_local < nx) {
	size_t iindex= (iix_local*nc + iy)*nckz; // iz=0
	double delta2= vol*fk[iindex][0];
	double delta_k_mag= 0.0;
	if(delta2 > 0.0)
	  delta_k_mag= sqrt(ampl*delta2);

	// delta(-k) = delta(k)^*
	fk[iindex][0]= delta_k_mag*cos(phase);
	fk[iindex][1]= -delta_k_mag*sin(phase);
      }
    }
  }
  
  gsl_rng_free(rng);

  //fprintf(stderr, "P_min= %e\n", P_min);
  //fprintf(stderr, "negative P(k): %zu\n", negative);  
}

void lognormal_compute_velogicy_grid(Grid const * const grid,
				     const int axis,
				     const double redshift,
				     const double omega_m,
				     Grid* const grid_v)
{
  assert(grid->mode == grid_mode_k);
  assert(0 <= axis && axis < 3);
  grid_v->mode= grid_mode_k;
  

  const int nc= grid->nc;
  const int nx= grid->local_nx;
  const int ix0= grid->local_x0;
  const double boxsize= grid->boxsize;
  
  const int nckz= nc/2 + 1;
  const double f= growth_f(1.0/(1.0 + redshift), omega_m);
  const double fac= f*boxsize/(2.0*M_PI);

  fftw_complex const * const d= grid->fk;
  fftw_complex * const v= grid_v->fk;


  double ik[3];
  
  for(int ix_local=0; ix_local<nx; ++ix_local) {
    int ix= ix0 + ix_local;
    ik[0]= ix < nc/2 ? ix : ix - nc;
    for(int iy=0; iy<nc; ++iy) {
      ik[1]= iy < nc/2 ? iy : iy - nc;
      
      int iz0 = ix == 0 && iy == 0; // skip kx=ky=kz=0
      for(int iz=iz0; iz<nckz; ++iz) {
	size_t index= (ix*static_cast<size_t>(nx) + iy)*nckz + iz;
	ik[2]= iz;

	double kfac= fac*(ik[axis]/(ik[0]*ik[0] + ik[1]*ik[1] + ik[2]*ik[2]));
	v[index][0]= kfac*d[index][1];
	v[index][1]= kfac*d[index][0]; // check negative sign
      }
    }
  }


}


void lognormal_generate_particles_periodic(const size_t np,
					   gsl_rng* rng,
					   Grid const * const grid_n,
					   Grid const * const grid_vx,
					   Grid const * const grid_vy,
					   Grid const * const grid_vz,
					   Particles& v)
{
  //
  // Add np particles (in average) to Particles
  //   grid_n: Grid of n(x), number of particles per (1/h Mpc)^3
  //   grid_vx: Grid of v_x(x) in 1/h Mpc (RSD displacement)
  //            This can be NULL (p.v[k] will be 0.0)
  //
  assert(grid_n->mode == grid_mode_x);
  if(grid_vx) {
    assert(grid_vx->mode == grid_mode_x);
    assert(grid_vy->mode == grid_mode_x);
    assert(grid_vz->mode == grid_mode_x);
  }
    
  const size_t nc= grid_n->nc;
  const size_t nx= grid_n->local_nx;
  const size_t ncz= 2*(nc/2 + 1);

  double const * const n= grid_n->fx;
  
  const double boxsize= grid_n->boxsize;
  v.boxsize= boxsize;

  const double num_bar= static_cast<double>(np)/(nc*nc*nc);
  const double dx= boxsize/nc;

  Particle p;
  p.v[0]= p.v[1]= p.v[2]= 0.0;
  //size_t np_added= 0;
  
  for(size_t ix_local=0; ix_local<nx; ++ix_local) {
    size_t ix= grid_n->local_x0 + ix_local;
    for(size_t iy=0; iy<nc; ++iy) {
      for(size_t iz=0; iz<nc; ++iz) {
	size_t index= (ix_local*nx + iy)*ncz + iz;

	// mean number of particles in the cell
	double num_grid= n[index]*num_bar;
	//num_total_mean += num_grid;
	//np_added++;

	// velocity
	if(grid_vx) {
	  p.v[0]= grid_vx->fx[index];
	  p.v[1]= grid_vy->fx[index];
	  p.v[2]= grid_vz->fx[index];	  
	}

	// poisson sampling of the number of particles in cell
	int num= gsl_ran_poisson(rng, num_grid);

	for(int i=0; i<num; ++i) {
	  p.x[0]= (ix + gsl_rng_uniform(rng))*dx;
	  p.x[1]= (iy + gsl_rng_uniform(rng))*dx;
	  p.x[2]= (iz + gsl_rng_uniform(rng))*dx;

	  v.push_back(p);
	}
      }
    }
  }

  //msg_printf(msg_debug, "%llu particles generated\n", np_added);
}

//
// Local functions
//
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
