#include "QMCFiniteSize/fsc_routines.h"

// =======================   basic I/O            =======================
xmlNodePtr find(const char* expression, xmlXPathContextPtr context)
{ // delegate to OhmmsXPathObject class
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

// =======================   structure manipulation   =======================
int get_grid_index1d(Ugrid grid, RealType x)
{
  return nearbyint((x-grid.start)/grid.delta);
}
int get_index3d_flat(Ugrid3D grid3d, int ix, int iy, int iz)
{
  return ix*grid3d.y.num*grid3d.z.num+iy*grid3d.z.num+iz;
}
