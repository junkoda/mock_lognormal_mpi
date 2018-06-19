#ifndef INPUT_POWER_H
#define INPUT_POWER_H 1

#include <vector>
//#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>

class InputPower {
 public:
  InputPower(const char filename[]);
  ~InputPower();
  double P(const double k) const;

 private:
  void load_txt(const char filename[]);
  void bcast();
  void init_interp();
  
  std::vector<double> log_k, log_P;
  gsl_interp* interp;
  gsl_interp_accel* acc;
  double k_min, k_max;
};

#endif
