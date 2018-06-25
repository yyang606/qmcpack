#ifndef QMCPLUSPLUS_NATURAL_SPLINE3D_H
#define QMCPLUSPLUS_NATURAL_SPLINE3D_H
#include "einspline/bspline_base.h"
#include "einspline/bspline_structs.h"
#include "einspline/bspline_create.h"
#include "einspline/bspline_eval_d.h"
namespace qmcplusplus
{
struct Ugrid3D
{
  Ugrid x, y, z;
};

class NaturalSpline3D
{
 public:
  NaturalSpline3D(const Ugrid3D grid3d, double* vals);
  ~NaturalSpline3D();
  BCtype_d natural_boundary();
  bool in_convex_hull(double x, double y, double z);

  // goal in life: evaluate Bspline at a point in 3D space
  double operator()(double x, double y, double z);

 private:
  const Ugrid3D grid3d_;
  UBspline_3d_d* spline3d_;  // major data structure to maintain
};
} // qmcplusplus
#endif
