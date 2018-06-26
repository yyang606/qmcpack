#include "QMCFiniteSize/fsc_routines.h"
#include "QMCFiniteSize/Spline3DFactory.h"
#include "QMCFiniteSize/SphericalAverage3D.h"
using namespace qmcplusplus;

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
  box.SetLRCutoffs();

  // step 2: obtain spline on regular grid (in reciprocal lattice units)
  vector<vector<RealType>> mat = loadtxt(fdat_name);
  Spline3DFactory spline_factory;
  string grid_name = "input_grid";
  Ugrid3D grid3d = spline_factory.create_ugrid3d(doc, grid_name);
  NaturalSpline3DInBox* boxspl3d = spline_factory.create_boxspl3d(
    grid3d, mat, box);

  // step 3: setup integral parameters
  int nrule=4, nk=64;
  RealType kmax = box.LR_kc;
  xmlNodePtr node;
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
  ofstream ofx, ofy, ofz;
  ofx.open("s100.dat");
  ofy.open("s110.dat");
  ofz.open("s111.dat");
  for (int ik=0; ik<nk; ik++)
  {
    RealType kmag = kmin+ik*dk;
    ofx << kmag << " " << (*boxspl3d)(kmag, 0, 0) << endl;
    ofy << kmag << " " << (*boxspl3d)(kmag/std::sqrt(2), kmag/std::sqrt(2), 0) << endl;
    ofz << kmag << " " << (*boxspl3d)(kmag/std::sqrt(3), kmag/std::sqrt(3), kmag/std::sqrt(3)) << endl;
  }
  ofx.close();
  ofy.close();
  ofz.close();

  OHMMS::Controller->finalize();
}
