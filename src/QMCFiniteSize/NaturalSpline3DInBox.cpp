#include "QMCFiniteSize/NaturalSpline3DInBox.h"

NaturalSpline3DInBox::NaturalSpline3DInBox(
  const NaturalSpline3D &spline3d,
  const Uniform3DGridLayout &box
) : spline3d_(spline3d), box_(box)
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
