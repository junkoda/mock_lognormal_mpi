#ifndef LOGNORMAL_H
#define LOGNORMAL_H 1

#include "grid.h"
#include "config.h"
#include "input_power.h"

Grid* lognormal_create_power_grid(InputPower const * const ps,
				  const int nc, const double boxsize,
				  const double pk_fac);

void lognormal_compute_gaussian_power(Grid* const grid);

void lognormal_create_gaussian_delta_k(const unsigned long seed,
				       Grid* const grid);

void lognormal_generate_particles_periodic(const size_t np,
					   gsl_rng* rng,
					   Grid const * const grid_n,
					   Grid const * const grid_vx,
					   Grid const * const grid_vy,
					   Grid const * const grid_vz,
					   Particles& v);

#endif
