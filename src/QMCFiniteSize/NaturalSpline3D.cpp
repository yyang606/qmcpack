#include "NaturalSpline3D.h"
namespace qmcplusplus
{
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
bool NaturalSpline3D::in_convex_hull(double x, double y, double z)
{
  bool inx, iny, inz;
  inx = iny = inz = true;
  if (x<grid3d_.x.start | x>grid3d_.x.end) inx = false;
  if (y<grid3d_.y.start | y>grid3d_.y.end) iny = false;
  if (z<grid3d_.z.start | z>grid3d_.z.end) inz = false;
  return (inx & iny & inz);
}
double NaturalSpline3D::operator()(double x, double y, double z)
{
  if (not in_convex_hull(x, y, z)) return NAN;
  double val;
  eval_UBspline_3d_d(spline3d_, x, y, z, &val);
  return val;
}
} // qmcplusplus
