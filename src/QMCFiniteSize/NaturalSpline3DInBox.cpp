#include "QMCFiniteSize/NaturalSpline3DInBox.h"
namespace qmcplusplus
{
NaturalSpline3DInBox::NaturalSpline3DInBox(
  const Ugrid3D grid3d,
  double* vals,
  const Uniform3DGridLayout &box
) : spline3d_(NaturalSpline3D(grid3d, vals)), box_(box)
{
}
NaturalSpline3DInBox::~NaturalSpline3DInBox()
{
}
double NaturalSpline3DInBox::operator()(double x, double y, double z)
{
  PosType kvec(x, y, z);
  PosType gvec = box_.k_unit(kvec);
  double myx = gvec[0];
  double myy = gvec[1];
  double myz = gvec[2];
  return spline3d_(myx, myy, myz);
}
} // qmcplusplus
