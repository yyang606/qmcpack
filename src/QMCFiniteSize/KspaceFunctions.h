#ifndef QMCPLUSPLUS_KSPACE_FUNCTIONS_H
#define QMCPLUSPLUS_KSPACE_FUNCTIONS_H
#define _USE_MATH_DEFINES
#include <cmath>
#include "QMCFiniteSize/BreakBase.h"
namespace qmcplusplus
{
class EslerCoul
{
 public:
  EslerCoul(RealType volume) : volume_(volume){};
  ~EslerCoul(){};
  RealType operator()(RealType k, RealType rc)
  {
    return 4.*M_PI/volume_/(k*k)*std::cos(k*rc);
  };
  void set_volume(RealType volume){ volume_=volume;};
 private:
  RealType volume_;
};
class NatoliCoul
{
 public:
  NatoliCoul(RealType volume) : volume_(volume){};
  ~NatoliCoul(){};
  RealType operator()(RealType k)
  {
    return 4.*M_PI/volume_/(k*k);
  };
  void set_volume(RealType volume){ volume_=volume;};
 private:
  RealType volume_;
};
} // qmcplusplus
#endif
