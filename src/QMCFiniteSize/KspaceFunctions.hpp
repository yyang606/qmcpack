#ifndef KSPACE_FUNCTIONS_H
#define KSPACE_FUNCTIONS_H
#define _USE_MATH_DEFINES
#include <cmath>

template <typename REAL>
class EslerCoul
{
 public:
  EslerCoul(REAL volume) : volume_(volume){};
  ~EslerCoul(){};
  double operator()(double k, double rc)
  {
    return 4.*M_PI/volume_*std::cos(k*rc)/(k*k);
  };
  void set_volume(REAL volume){ volume_=volume;};
 private:
  REAL volume_;
};
#endif
