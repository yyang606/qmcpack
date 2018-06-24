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
Uniform3DGridLayout create_box(xmlXPathContextPtr doc)
{ // delegate to LatticeParser clasee
  xmlNodePtr sc_node = find("//simulationcell", doc);

  Uniform3DGridLayout box;
  LatticeParser parser(box);
  parser.put(sc_node);
  return box;
}

// =======================   structure manipulation   =======================
// -----------------------   spline             -----------------------
vector<RealType> spherical_integral(
  NaturalSpline3DInBox boxspl3d,
  vector<RealType> kmags,
  int nrule
)
{
  // integral results
  vector<RealType> intvals(kmags.size(), 0.0);

  // initialize quadrature points and weights
  Quadrature3D<RealType> qrule(nrule);

  for (int ik=0; ik<kmags.size(); ik++)
  {
    RealType kval = kmags[ik];
    for (int i=0; i<qrule.xyz_m.size(); i++)
    {
      double kx = kval*qrule.xyz_m[i][0];
      double ky = kval*qrule.xyz_m[i][1];
      double kz = kval*qrule.xyz_m[i][2];
      intvals[ik] += boxspl3d(kx, ky, kz)*qrule.weight_m[i];
    }
  }
  return intvals;
}
