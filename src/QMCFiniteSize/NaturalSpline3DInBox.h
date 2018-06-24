#ifndef QMCPLUSPLUS_NATURAL_SPLINE3D_IN_BOX_H
#define QMCPLUSPLUS_NATURAL_SPLINE3D_IN_BOX_H
#include "QMCFiniteSize/NaturalSpline3D.h"
#include "Lattice/Uniform3DGridLayout.h"
#include "Configuration.h"
namespace qmcplusplus
{
typedef QMCTraits::PosType PosType;
class NaturalSpline3DInBox
{
 public:
  NaturalSpline3DInBox(
    const Ugrid3D grid3d,
    double* vals,
    const Uniform3DGridLayout &box
  );
  ~NaturalSpline3DInBox();
  double operator()(double x, double y, double z);

 private:
  NaturalSpline3D spline3d_;
  const Uniform3DGridLayout box_;
};
} // qmcplusplus
#endif
