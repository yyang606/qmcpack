#ifndef FSC_COMMON_H
#define FSC_COMMON_H
#include <fstream>
#include <sstream>
#include <vector>
#include "OhmmsData/Libxml2Doc.h"
#include "Lattice/Uniform3DGridLayout.h"
#include "ParticleIO/ParticleLayoutIO.h"

using namespace qmcplusplus;
using namespace std;

typedef QMCTraits::RealType  RealType;

// =======================   basic I/O            =======================

// Finds the first match for the given xpath expression in given context.
//
// Example:
//   xmlNodePtr sc_node = find("//simulationcell", doc);
//   xmlNodePtr bc_node = find('//parameter[@name="bconds"]', sc_node);
xmlNodePtr find(const char* expression, xmlXPathContextPtr context);

// Reads a table of numbers into a matrix. Essentially np.loadtxt in Python
//  for a 2D table of real numbers.
//
// Example:
//   vector<vector<RealType>> mat = loadtxt("matrix.dat");
vector<vector<RealType>> loadtxt(const string fname);

// =======================   structure creation   =======================
// -----------------------   box                  -----------------------
Uniform3DGridLayout create_box(xmlNodePtr sc_node);
#endif
