#include "QMCFiniteSize/fsc_routines.h"

// =======================   basic I/O            =======================
xmlNodePtr find(const char* expression, xmlXPathContextPtr context)
{ // delegate to OhmmsXPathObject class
  OhmmsXPathObject xpath(expression, context);
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
Uniform3DGridLayout create_box(xmlNodePtr sc_node)
{ // delegate to LatticeParser clasee
  Uniform3DGridLayout box;
  LatticeParser parser(box);
  parser.put(sc_node);
  return box;
}
