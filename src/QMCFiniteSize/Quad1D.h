#ifndef QUAD1D_H
#define QUAD1D_H
#include <vector>
using namespace std;
class Quad1D
{
 public:
  vector<double> x, w;  // abscissas and weights
  Quad1D(double a, double b, int nk);
  double integrate(vector<double> y);
};
#endif
