#ifndef SPLINE1D_H
#define SPLINE1D_H
#include <vector>
#include "einspline/nubspline_create.h"
#include "einspline/nubspline_eval_d.h"
class Spline1D
{
 public:
  Spline1D(std::vector<double> pts, std::vector<double> vals);
  Spline1D(std::vector<double> pts, std::vector<double> vals, BCtype_d bc);
  double operator()(double x);
 private:
  NUgrid *grid1d_;
  NUBspline_1d_d *spline1d_;
};
#endif
