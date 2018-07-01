#ifndef LOGNORMAL_H
#define LOGNORMAL_H 1

#include "grid.h"
#include "input_power.h"

Grid* lognormal_create_power_grid(InputPower const * const ps,
				  const int nc, const double boxsize,
				  const double pk_fac);

void lognormal_compute_gaussian_power(Grid* const grid);

void lognormal_create_gaussian_delta_k(const unsigned long seed,
				       Grid* const grid);

#endif
