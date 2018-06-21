#include "QMCFiniteSize/fsc_routines.h"

int main(int argc, char **argv)
{
  OHMMS::Controller->initialize(argc, argv);
  if (argc<3) APP_ABORT("usage: " << argv[0] << " axes.xml sofk.dat");

  string fxml_name = argv[1];
  string fdat_name = argv[2];
 
  // parse input xml
  Libxml2Document fxml;
  fxml.parse(fxml_name);
  xmlXPathContextPtr doc = fxml.getXPathContext();

  // step 1: construct lattice -> enable box.k_unit, box.k_cart
  Uniform3DGridLayout box = create_box(doc);

  // step 2: load data on simulation kgrid
  vector<vector<RealType>> mat = loadtxt(fdat_name);
  int nk = mat.size();

  // step 3: initialize regular grid
  string grid_name = "input_grid";
  Ugrid3D grid3d = create_ugrid3d(doc, grid_name);
  int nval = grid3d.x.num*grid3d.y.num*grid3d.z.num;
  if (nval!=nk)
  {
    app_log() << "expected " << nval << " found " << nk << endl;
    APP_ABORT("grid mismatch");
  }

  // step 4: transfer data to regular grid
  RealType* vals = new RealType[nval];
  app_log() << "transfer data to regular grid" << endl;
  for (int ik=0; ik<nk; ik++)
  {
    PosType kvec(mat[ik][0], mat[ik][1], mat[ik][2]);
    PosType gvec = box.k_unit(kvec);
    int ix = get_grid_index1d(grid3d.x, gvec[0]);
    int iy = get_grid_index1d(grid3d.y, gvec[1]);
    int iz = get_grid_index1d(grid3d.z, gvec[2]);
    int idx = get_index3d_flat(grid3d, ix, iy, iz);
    vals[idx] = mat[ik][3];
  }
  NaturalSpline3D spline3d = NaturalSpline3D(grid3d, vals);
  delete vals;
  NaturalSpline3DInBox boxspl3d(spline3d, box);

  // step 5: dump regular grid and spline for debugging
  app_log() << "dump regular grid and spline" << endl;
  ofstream ofs("grid.dat", ofstream::out);
  for (int ix=0; ix<grid3d.x.num; ix++)
  {
    for (int iy=0; iy<grid3d.y.num; iy++)
    {
      for (int iz=0; iz<grid3d.z.num; iz++)
      {
        int idx = get_index3d_flat(grid3d, ix, iy, iz);
        RealType gx = get_grid_point1d(grid3d.x, ix);
        RealType gy = get_grid_point1d(grid3d.y, iy);
        RealType gz = get_grid_point1d(grid3d.z, iz);
        //RealType val = spline3d(gx, gy, gz);
        PosType gvec(gx, gy, gz);
        PosType kvec = box.k_cart(gvec);
        RealType kx = kvec[0];
        RealType ky = kvec[1];
        RealType kz = kvec[2];
        RealType val = boxspl3d(kx, ky, kz);
        ofs << gx << " " << gy << " " << gz << " "
            << val << endl;
      }
    }
  }
  ofs.close();

  int nrule = 4;
  nk = 64;
  RealType kmax = 2.5;
  RealType dk = kmax/nk;
  vector<RealType> kmags(nk);
  for (int ik=0; ik<nk; ik+=1)
  {
    kmags[ik] = ik*dk;
  }
  vector<RealType> intvals = spherical_integral(boxspl3d, kmags, nrule);
  ofs.open("avesk.dat", ofstream::out);
  for (int ik=0; ik<nk; ik+=1)
  {
    ofs << kmags[ik] << " " << intvals[ik] << endl;
  }
  ofs.close();

}
