#include <fstream>
#include <sstream>
#include "QMCFiniteSize/fsc_routines.h"
#include "einspline/bspline_create.h"

// =======================   basic I/O            =======================

xmlNodePtr find(const char* expression, xmlXPathContextPtr context)
{ // find the first node matching xpath expression in the given context
  OhmmsXPathObject xpath(expression, context);
  if (xpath.size() != 1)
  {
    APP_ABORT("expected 1 " << expression << " found " << xpath.size());
  }
  return xpath[0];
}


vector<vector<RealType>> loadtxt(const string fname)
{ // stack overflow "How to read in a data file of unknown dimensions in C/C++"

  vector<vector<RealType>> mat;
  ifstream fdat;
  fdat.open(fname.c_str(), ios_base::in);

  double buf;
  string line;

  while (!fdat.eof())
  {
    getline(fdat, line);

    // skip comment or empty line
    if (line[0] == '#' || line.empty()) continue;

    // read a row of numbers
    vector<RealType> row;
    stringstream ss(line, ios_base::in);
    while (ss>>buf) row.push_back(buf);

    // add a row to matrix
    mat.push_back(row);
  }

  fdat.close();
  return mat;
}

// =======================   structure creation   =======================

// -----------------------   box                  -----------------------

Uniform3DGridLayout create_box(xmlXPathContextPtr doc)
{
  xmlNodePtr sc_node = find("//simulationcell", doc);

  Uniform3DGridLayout box;
  LatticeParser parser(box);
  parser.put(sc_node);
  return box;
}


BCtype_d natural_boundary()
{
  BCtype_d bc;
  bc.lCode = NATURAL;
  bc.rCode = NATURAL;
  bc.lVal = 1.0;
  bc.rVal = 1.0;
  return bc;
}

// -----------------------   spline               -----------------------

UBspline_3d_d* create_natural_spline3d(
  GridType gridx,
  GridType gridy,
  GridType gridz,
  RealType* data)
{ // spline 3D scalar field on regular grid
  Ugrid esgridx, esgridy, esgridz;
  esgridx = gridx.einspline_grid();
  esgridy = gridy.einspline_grid();
  esgridz = gridz.einspline_grid();
  BCtype_d bcx, bcy, bcz;
  bcx = natural_boundary();
  bcy = natural_boundary();
  bcz = natural_boundary();
  UBspline_3d_d* spline3d = create_UBspline_3d_d(
    esgridx, esgridy, esgridz
  , bcx, bcy, bcz, data
  );
  return spline3d;
}
