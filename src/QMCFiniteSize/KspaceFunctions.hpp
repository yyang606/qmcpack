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
  REAL operator()(REAL k, REAL rc)
  {
    return 4.*M_PI/volume_/(k*k)*std::cos(k*rc);
  };
  void set_volume(REAL volume){ volume_=volume;};
 private:
  REAL volume_;
};
template <typename REAL>
class NatoliCoul
{
 public:
  NatoliCoul(REAL volume) : volume_(volume){};
  ~NatoliCoul(){};
  REAL operator()(REAL k)
  {
    return 4.*M_PI/volume_/(k*k);
  };
  void set_volume(REAL volume){ volume_=volume;};
 private:
  REAL volume_;
};
#endif
