#include <cstdio>
#include <cmath>
#include <cassert>

#include "comm.h"
#include "msg.h"
#include "error.h"
#include "input_power.h"


InputPower::InputPower(const char filename[])
{
  if(comm_this_node() == 0)
    load_txt(filename);

  bcast();
  init_interp();
}

InputPower::~InputPower()
{
  gsl_interp_accel_free(acc);
  gsl_interp_free(interp);
}

void InputPower::load_txt(const char filename[])
{
  //
  // Load k P pairs from ascii file and store in log_k, log_P
  //
  FILE* fp= fopen(filename, "r");
  if(fp == 0) {
    throw FileNotFoundError();
  }

  double k, P;

  char line[1024];
  
  while(fgets(line, 1024, fp)) {
    if(line[0] == '#')
      continue;

    int ret = sscanf(line, "%lg %lg", &k, &P);

    if(ret != 2) {
      msg_abort("Error: unable to parse line: %s", line);
    }

    log_k.push_back(log(k));
    log_P.push_back(log(P));
  }

  int ret= fclose(fp); assert(ret == 0);
  assert(log_k.size() == log_P.size());

  if(log_k.size() == 0) {
    msg_abort("No lines read from: %s", filename);
  }
}

void InputPower::bcast()
{
  // Share log_k, log_P with other MPI nodes
  int n= static_cast<int>(log_k.size());

  comm_bcast(n);

  if(comm_this_node() != 0) {
    log_k.resize(n, 0.0);
    log_P.resize(n, 0.0);
  }
  
  MPI_Bcast(log_k.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(log_P.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void InputPower::init_interp()
{
  interp= gsl_interp_alloc(gsl_interp_cspline, log_k.size());
  acc= gsl_interp_accel_alloc();

  const size_t n_required= gsl_interp_min_size(interp);
  if(log_k.size() < n_required) {
    msg_abort("Error: Not enough power spectrum data points for cubic spline; %lu data points < %lu required\n", log_k.size(), n_required);
  }

  gsl_interp_init(interp, log_k.data(), log_P.data(), log_k.size());

  k_min= exp(log_k[0]);
  assert(log_k.size() >= 1);
  k_max= exp(log_k[log_k.size()-1]);
}

double InputPower::P(const double k) const
{
  if(k < k_min || k > k_max) {
    msg_abort("k is beyond the interpolation range: %e %e %e\n",
	      k, k_min, k_max);
  }
  
  double log_Pk=
    gsl_interp_eval(interp, log_k.data(), log_P.data(), log(k), acc);
  
  return exp(log_Pk);
}
