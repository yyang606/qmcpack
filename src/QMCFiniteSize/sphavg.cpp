#include "QMCFiniteSize/fsc_routines.h"
#include "QMCFiniteSize/Spline3DFactory.h"
#include "QMCFiniteSize/SphericalAverage3D.h"
using namespace qmcplusplus;

vector<RealType> get_hessian_at_zero(NaturalSpline3DInBox* boxspl3d, RealType h)
{
  RealType h2 = h*h;
  RealType x, y, z;
  x=0; y=0; z=0;

  // evaluate spline at stencils
  RealType f, fx_p, fx_m, fy_p, fy_m, fz_p, fz_m;
  RealType fxy_p, fxy_m, fxz_p, fxz_m, fyz_p, fyz_m;
  f    = (*boxspl3d)(x, y, z);
  fx_p = (*boxspl3d)(x+h, y, z);
  fx_m = (*boxspl3d)(x-h, y, z);
  fy_p = (*boxspl3d)(x, y+h, z);
  fy_m = (*boxspl3d)(x, y-h, z);
  fz_p = (*boxspl3d)(x, y, z+h);
  fz_m = (*boxspl3d)(x, y, z-h);
  fxy_p = (*boxspl3d)(x+h, y+h, z);
  fxy_m = (*boxspl3d)(x-h, y-h, z);
  fxz_p = (*boxspl3d)(x+h, y, z+h);
  fxz_m = (*boxspl3d)(x-h, y, z-h);
  fyz_p = (*boxspl3d)(x, y+h, z+h);
  fyz_m = (*boxspl3d)(x, y-h, z-h);

  // build hessian matrix
  RealType hxx, hyy, hzz, hxy, hxz, hyz;
  hxx = (fx_p+fx_m-2*f)/h2;
  hyy = (fy_p+fy_m-2*f)/h2;
  hzz = (fz_p+fz_m-2*f)/h2;
  hxy = (fxy_p+fxy_m+2*f-fx_p-fx_m-fy_p-fy_m)/(2*h2);
  hxz = (fxz_p+fxz_m+2*f-fx_p-fx_m-fz_p-fz_m)/(2*h2);
  hyz = (fyz_p+fyz_m+2*f-fy_p-fy_m-fz_p-fz_m)/(2*h2);

  vector<RealType> hessian(6);
  hessian[0] = hxx;
  hessian[1] = hyy;
  hessian[2] = hzz;
  hessian[3] = hxy;
  hessian[4] = hxz;
  hessian[5] = hyz;
  return hessian;
}

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
  xmlNodePtr sc_node = find("//simulationcell", doc);
  Uniform3DGridLayout box = create_box(sc_node);
  box.SetLRCutoffs();

  // step 2: obtain spline on regular grid (in reciprocal lattice units)
  vector<vector<RealType>> mat = loadtxt(fdat_name);
  Spline3DFactory spline_factory;
  string grid_name = "input_grid";
  Ugrid3D grid3d = spline_factory.create_ugrid3d(doc, grid_name);
  NaturalSpline3DInBox* boxspl3d = spline_factory.create_boxspl3d(
    grid3d, mat, box);

  xmlNodePtr node;
  RealType h = 0.0001;
  node = find("//hessian/h", doc);
  if (node) putContent(h, node);

  // !!!! get hessian matrix
  vector<RealType> hess = get_hessian_at_zero(boxspl3d, h);
  ofstream ofh;
  ofh.open("hess.dat");
  for (int i=0;i<hess.size();i++)
  {
    ofh << hess[i] << endl;
  }
  ofh.close();

  // step 3: setup integral parameters
  int nrule=4, nk=64;
  RealType kmax = box.LR_kc;
  node = find("//output_grid/nk", doc);
  if (node) putContent(nk, node);
  node = find("//output_grid/kmax", doc);
  if (node) putContent(kmax, node);
  node = find("//spherical_average/nrule", doc);
  if (node) putContent(nrule, node);
  app_log() << " nrule = " << nrule << endl;
  SphericalAverage3D sphavg(nrule);

  // output useful stuff
  RealType dk = kmax/(nk-1);
  ofstream ofs;
  ofs.open("avesk.dat");
  for (int ik=0; ik<nk; ik+=1)
  {
    RealType kmag = ik*dk;
    RealType sk = sphavg(*boxspl3d, kmag);
    ofs << kmag << " " << sk << endl;
  }
  ofs.close();

  RealType kmin = 0.0;
  ofstream ofx, ofy, ofz, of1, of2, of3;
  of1.open("s100.dat");
  of2.open("s110.dat");
  of3.open("s111.dat");
  ofx.open("sx.dat");
  ofy.open("sy.dat");
  ofz.open("sz.dat");
  for (int ik=0; ik<nk; ik++)
  {
    RealType kmag = kmin+ik*dk;
    of1 << kmag << " " << (*boxspl3d)(kmag, 0, 0) << endl;
    of2 << kmag << " " << (*boxspl3d)(kmag/std::sqrt(2), kmag/std::sqrt(2), 0) << endl;
    of3 << kmag << " " << (*boxspl3d)(kmag/std::sqrt(3), kmag/std::sqrt(3), kmag/std::sqrt(3)) << endl;
    ofx << kmag << " " << (*boxspl3d)(kmag, 0, 0) << endl;
    ofy << kmag << " " << (*boxspl3d)(0, kmag, 0) << endl;
    ofz << kmag << " " << (*boxspl3d)(0, 0, kmag) << endl;
  }
  of1.close();
  of2.close();
  of3.close();
  ofx.close();
  ofy.close();
  ofz.close();

  OHMMS::Controller->finalize();
}
