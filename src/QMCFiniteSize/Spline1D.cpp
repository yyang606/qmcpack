#include "Spline1D.h"
using namespace std;
Spline1D::Spline1D(vector<double> pts, vector<double> vals)
{
  grid1d_ = create_general_grid(pts.data(), pts.size());
  BCtype_d bc;
  bc.lVal = 0.0;
  bc.rVal = 0.0;
  bc.lCode=DERIV1;
  bc.rCode=DERIV1;
  spline1d_ = create_NUBspline_1d_d(grid1d_, bc, vals.data());
}
double Spline1D::operator()(double x)
{
  double val;
  eval_NUBspline_1d_d(spline1d_, x, &val);
  return val;
}
