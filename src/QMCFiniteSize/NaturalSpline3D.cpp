#include "NaturalSpline3D.h"

NaturalSpline3D::NaturalSpline3D(const Ugrid3D grid3d, double* vals)
  : grid3d_(grid3d)
{
  BCtype_d bcx, bcy, bcz;
  bcx = natural_boundary();
  bcy = natural_boundary();
  bcz = natural_boundary();
  spline3d_ = create_UBspline_3d_d(
    grid3d.x, grid3d.y, grid3d.z
  , bcx, bcy, bcz, vals
  );
}

NaturalSpline3D::~NaturalSpline3D()
{
  delete spline3d_;
}

BCtype_d NaturalSpline3D::natural_boundary()
{
  BCtype_d bc;
  bc.lCode = NATURAL;
  bc.rCode = NATURAL;
  bc.lVal = 1.0;
  bc.rVal = 1.0;
  return bc;
}

double NaturalSpline3D::operator()(double x, double y, double z)
{
  double val;
  eval_UBspline_3d_d(spline3d_, x, y, z, &val);
  return val;
}
