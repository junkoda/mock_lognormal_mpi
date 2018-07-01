#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>


//static double D0= 0.0;  // normalization factor

static double growth_int(double a, void* param);
static double growth(const double a, const double omega_m);

double growth_D(const double a, const double omega_m)
{
  // Growth factor D(a)
  return growth(a, omega_m)/growth(1.0, omega_m);
}

double growth_f(const double a, const double omega_m)
{
  // Growth rate f=dlnD/dlna
  const double D0= growth(1.0, omega_m);
  const double hubble_a= sqrt(omega_m/(a*a*a) + 1.0 - omega_m);
  const double d= growth(a, omega_m)/D0;

  const double f_ex= 1.0/(d*D0*a*a*hubble_a*hubble_a)
                     - 1.5*omega_m/(hubble_a*hubble_a*a*a*a);

  return f_ex;
}



double growth_int(double a, void* param)
{
  const double om= *(double*)param;
  return pow(a/(om + (1.0 - om)*a*a*a), 1.5);
}

double growth(const double a, const double omega_m)
{
  // Compute integral of D_tilda(a) -- unnormalised growth factor
  // growth_factor = D_tilda(a)/D_tilda(1.0)

  const double hubble_a= sqrt(omega_m/(a*a*a) + 1.0 - omega_m);

  const int worksize= 1000;
  double result, abserr;

  gsl_integration_workspace *workspace=
    gsl_integration_workspace_alloc(worksize);

  gsl_function F;
  F.function = &growth_int;
  F.params = (void*) &omega_m;

  gsl_integration_qag(&F, 0, a, 0, 1.0e-8, worksize, GSL_INTEG_GAUSS41,
		      workspace, &result, &abserr);

  gsl_integration_workspace_free(workspace);

  return hubble_a * result;
}
