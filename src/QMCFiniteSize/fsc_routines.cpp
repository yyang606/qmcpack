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
// -----------------------   grid               -----------------------
Ugrid create_ugrid1d(xmlNodePtr node)
{
  vector<double> min_max_num;
  putContent(min_max_num, node);
  Ugrid grid;
  grid.start = min_max_num[0];
  grid.end   = min_max_num[1];
  grid.num   = (int) nearbyint(min_max_num[2]);
  grid.delta = (grid.end-grid.start)/(grid.num-1);
  grid.delta_inv = 1./grid.delta;
  return grid;
}
Ugrid3D create_ugrid3d(xmlXPathContextPtr doc, string name)
{ // use create_ugrid1d to make 3 1D grids, stuff'em in a 3D grid
  string expr = "//ugrid3d[\"" + name + "\"]";
  xmlNodePtr ug_node = find(expr.c_str(), doc);

  // initialize
  Ugrid3D grid3d;
  vector<bool> created(3, false);
  xmlNodePtr node=ug_node->children;
  while (node != NULL)
  {
    string name((const char*)node->name);
    if (name == "x")
    {
      grid3d.x = create_ugrid1d(node);
      created[0] = true;
    }
    if (name == "y")
    {
      grid3d.y = create_ugrid1d(node);
      created[1] = true;
    }
    if (name == "z")
    {
      grid3d.z = create_ugrid1d(node);
      created[2] = true;
    }
    node = node->next;
  }

  // check
  for (int idim=0; idim<3; idim++)
  {
    if (not created[idim])
    {
      APP_ABORT("create_ugrid3d failed");
    }
  }
  
  return grid3d;
}
// -----------------------   spline             -----------------------
NaturalSpline3DInBox create_boxspl3d(
  Ugrid3D grid3d, vector<vector<RealType>> mat, Uniform3DGridLayout box
)
{
  int nk = mat.size();
  int nval = grid3d.x.num*grid3d.y.num*grid3d.z.num;
  if (nval!=nk) APP_ABORT("grid mismatch: input " << nval << " data " << nk);

  // transfer data to regular grid
  vector<RealType> vals(nval);
  vector<bool> filled(nval, false);
  for (int ik=0; ik<nk; ik++)
  {
    PosType kvec(mat[ik][0], mat[ik][1], mat[ik][2]);
    PosType gvec = box.k_unit(kvec);
    int ix = get_grid_index1d(grid3d.x, gvec[0]);
    int iy = get_grid_index1d(grid3d.y, gvec[1]);
    int iz = get_grid_index1d(grid3d.z, gvec[2]);
    int idx = get_index3d_flat(grid3d, ix, iy, iz);
    vals[idx] = mat[ik][3];
    filled[idx] = true;
  }

  // verify all grid points are filled
  for (int ik=0; ik<nk; ik++)
  {
    if (not filled[ik]) APP_ABORT(ik << " is missing");
  }

  NaturalSpline3DInBox boxspl3d(grid3d, vals.data(), box);
  return boxspl3d;
}

// =======================   structure manipulation   =======================
// -----------------------   grid               -----------------------
RealType get_grid_point1d(Ugrid grid, int ix)
{
  return grid.start+ix*grid.delta;
}
int get_grid_index1d(Ugrid grid, RealType x)
{
  return nearbyint((x-grid.start)/grid.delta);
}
int get_index3d_flat(Ugrid3D grid3d, int ix, int iy, int iz)
{
  return ix*grid3d.y.num*grid3d.z.num+iy*grid3d.z.num+iz;
}
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
