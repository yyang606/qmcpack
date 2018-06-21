#include "QMCFiniteSize/NaturalSpline3D.h"
#include "Lattice/Uniform3DGridLayout.h"
#include "Configuration.h"

using namespace qmcplusplus;

typedef QMCTraits::PosType PosType;

class NaturalSpline3DInBox
{
 public:
  NaturalSpline3DInBox(
    const NaturalSpline3D &spline3d,
    const Uniform3DGridLayout &box
  );
  ~NaturalSpline3DInBox();
  double operator()(double x, double y, double z);
 private:
  NaturalSpline3D spline3d_;
  const Uniform3DGridLayout box_;
};
