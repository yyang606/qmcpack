#include "Spline1D.h"
using namespace std;
Spline1D::Spline1D(vector<double> pts, vector<double> vals)
{
  BCtype_d bc;
  bc.lVal = 0.0;
  bc.rVal = 0.0;
  bc.lCode=DERIV1;
  bc.rCode=DERIV1;
  Spline1D(pts, vals, bc);
}
Spline1D::Spline1D(vector<double> pts, vector<double> vals, BCtype_d bc)
{
  grid1d_ = create_general_grid(pts.data(), pts.size());
  spline1d_ = create_NUBspline_1d_d(grid1d_, bc, vals.data());
}
double Spline1D::operator()(double x)
{
  double val;
  eval_NUBspline_1d_d(spline1d_, x, &val);
  return val;
}
