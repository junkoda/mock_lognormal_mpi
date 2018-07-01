#ifndef CONFIG_H
#define CONFIG_H 1

#include <vector>

struct Particle {
  double x[3], v[3];
};

class Particles : public std::vector<Particle> {
 public:
  double boxsize;
};

#endif
