#ifndef QMCPLUSPLUS_SPHERICAL_AVERAGE3D_H
#define QMCPLUSPLUS_SPHERICAL_AVERAGE3D_H
#include "QMCFiniteSize/NaturalSpline3DInBox.h"
#include "Numerics/Quadrature.h"
#include "Configuration.h"
#include <functional>
namespace qmcplusplus
{
typedef QMCTraits::RealType RealType;
// Spherically average a 3D scalar field at given radius
class SphericalAverage3D
{
 public:
  SphericalAverage3D(int nrule);
  RealType operator()(
    std::function< RealType(RealType, RealType, RealType) > func,
    RealType r
  );
  void report(std::ostream &os);
 private:
  Quadrature3D<RealType> qrule_;
};
} // qmcplusplus
#endif
