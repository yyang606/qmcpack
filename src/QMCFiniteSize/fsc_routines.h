#include <vector>
#include "OhmmsData/Libxml2Doc.h"
#include "Numerics/OneDimGridBase.h"
#include "Lattice/Uniform3DGridLayout.h"
#include "ParticleIO/ParticleLayoutIO.h"
#include "einspline/bspline_base.h"
#include "einspline/bspline_structs.h"

using namespace qmcplusplus;
using namespace std;

typedef               double  RealType;
typedef LinearGrid<RealType>  GridType;

// =======================   basic I/O            =======================

xmlNodePtr find(const char* expression, xmlXPathContextPtr context);

vector<vector<RealType>> loadtxt(const string fname);

// =======================   structure creation   =======================

// -----------------------   box                  -----------------------
Uniform3DGridLayout create_box(xmlXPathContextPtr doc);

// -----------------------   spline               -----------------------

BCtype_d natural_boundary();

UBspline_3d_d* create_natural_spline3d(
  GridType gridx,
  GridType gridy,
  GridType gridz,
  RealType* data);
