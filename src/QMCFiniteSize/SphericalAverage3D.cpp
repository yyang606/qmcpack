#include "QMCFiniteSize/SphericalAverage3D.h"
using namespace std;
namespace qmcplusplus
{
SphericalAverage3D::SphericalAverage3D(int nrule) : qrule_(nrule){}
RealType SphericalAverage3D::operator()(
  function< RealType(RealType, RealType, RealType) > func3d,
  RealType r
)
{
  RealType intval = 0;
  for (int i=0; i<qrule_.xyz_m.size(); i++)
  {
    double x = r*qrule_.xyz_m[i][0];
    double y = r*qrule_.xyz_m[i][1];
    double z = r*qrule_.xyz_m[i][2];
    intval += func3d(x, y, z)*qrule_.weight_m[i];
  }
  return intval;
}
void SphericalAverage3D::report(ostream &os)
{
  for (int i=0; i<qrule_.xyz_m.size(); i++)
  {
    os << qrule_.xyz_m[i] << " " << qrule_.weight_m[i] << endl;
  }
}
} // qmcplusplus
