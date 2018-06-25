#include "Quad1D.h"
#include "quadrule.h"
Quad1D::Quad1D(double a, double b, int nk)
{
  x.resize(nk);
  w.resize(nk);
  // nk, kind (1: Gauss-Legendre), alpha, beta, a, b
  cgqf (nk, 1, 0.0, 0.0, a, b,
    -1, x.data(), w.data());
  // lo (<0: print, 0: compute, 1: print & check), t[], wgt[]
}
double Quad1D::integrate(vector<double> y)
{
  double sum=0;
  for (int ik=0; ik<x.size(); ik++)
  {
    sum += y[ik]*w[ik];
  }
  return sum;
}
