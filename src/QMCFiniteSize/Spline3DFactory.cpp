#include "QMCFiniteSize/Spline3DFactory.h"
using namespace std;
namespace qmcplusplus
{
NaturalSpline3DInBox* Spline3DFactory::create_boxspl3d(
  Ugrid3D grid3d,
  vector<vector<RealType>> mat,
  Uniform3DGridLayout box
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
  
  NaturalSpline3DInBox *boxspl3d = new 
    NaturalSpline3DInBox(grid3d, vals.data(), box);
  return boxspl3d;
}
Ugrid Spline3DFactory::create_ugrid1d(xmlNodePtr node)
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
Ugrid3D Spline3DFactory::create_ugrid3d(xmlXPathContextPtr doc, string name)
{ // use create_ugrid1d to make 3 1D grids, stuff'em in a 3D grid
  string expr = "//ugrid3d[@name=\"" + name + "\"]";
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
RealType Spline3DFactory::get_grid_point1d(Ugrid grid, int ix)
{
  return grid.start+ix*grid.delta;
}
int Spline3DFactory::get_grid_index1d(Ugrid grid, RealType x)
{
  return nearbyint((x-grid.start)/grid.delta);
}
int Spline3DFactory::get_index3d_flat(Ugrid3D grid3d, int ix, int iy, int iz)
{
  return ix*grid3d.y.num*grid3d.z.num+iy*grid3d.z.num+iz;
}
} // qmcplusplus
