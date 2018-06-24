#ifndef QMCPLUSPLUS_SPLINE3D_FACTORY_H
#define QMCPLUSPLUS_SPLINE3D_FACTORY_H
#include "QMCFiniteSize/fsc_routines.h"
#include "QMCFiniteSize/NaturalSpline3DInBox.h"
#include "Lattice/Uniform3DGridLayout.h"
#include "Configuration.h"
namespace qmcplusplus
{
typedef QMCTraits::RealType RealType;
typedef QMCTraits::PosType   PosType;
class Spline3DFactory
{
 public:
  NaturalSpline3DInBox* create_boxspl3d(
    Ugrid3D grid3d,
    std::vector<std::vector<RealType>> mat,
    Uniform3DGridLayout box
  );
  Ugrid   create_ugrid1d(xmlNodePtr node);
  Ugrid3D create_ugrid3d(xmlXPathContextPtr doc, std::string name);

 private:
  RealType get_grid_point1d(Ugrid grid, int ix);
  int get_grid_index1d(Ugrid grid, RealType x);
  int get_index3d_flat(Ugrid3D grid3d, int ix, int iy, int iz);
};
} // qmcplusplus
#endif
