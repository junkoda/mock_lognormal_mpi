#ifndef MASS_ASSIGNMENT_H
#define MASS_ASSIGNMENT_H 1

#include "grid.h"
#include "config.h"

void mass_assignment(const Particles& v,
		     const int mas,
		     Grid* const grid,
		     Grid* const grid_shifted);

#endif
