#include <fstream>
#include <sstream>
#include <vector>
#include "OhmmsData/Libxml2Doc.h"
#include "Lattice/Uniform3DGridLayout.h"
#include "ParticleIO/ParticleLayoutIO.h"
#include "QMCFiniteSize/NaturalSpline3D.h"

using namespace qmcplusplus;
using namespace std;

typedef QMCTraits::RealType  RealType;
typedef QMCTraits::PosType    PosType;

// =======================   basic I/O            =======================

// Finds the only match for the given xpath expression in given context.
//
// Example:
//   xmlNodePtr sc_node = find("//simulationcell", doc);
//   xmlNodePtr bc_node = find('//parameter[@name="bconds"]', sc_node);
xmlNodePtr find(const char* expression, xmlXPathContextPtr context);

// Reads a table of numbers into a matrix. Essentially np.loadtxt in Python
//  for a 2D table of floats.
//
// Example:
//   vector<vector<RealType>> mat = loadtxt("matrix.dat");
vector<vector<RealType>> loadtxt(const string fname);

// =======================   structure creation   =======================

// -----------------------   box                  -----------------------
Uniform3DGridLayout create_box(xmlXPathContextPtr doc);

// -----------------------   grid                 -----------------------
Ugrid   create_ugrid1d(xmlNodePtr ug_node);
Ugrid3D create_ugrid3d(xmlXPathContextPtr doc, string name);
int get_grid_index1d(Ugrid grid, RealType x);
int get_index3d_flat(Ugrid3D grid3d, int ix, int iy, int iz);
