#include "QMCFiniteSize/lr_routines.h"

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

  // step 2: obtain spline on regular grid (in reciprocal lattice units)
  vector<vector<RealType>> mat = loadtxt(fdat_name);
  string grid_name = "input_grid";
  Ugrid3D grid3d = create_ugrid3d(doc, grid_name);
  NaturalSpline3DInBox boxspl3d = create_boxspl3d(grid3d, mat, box);

  // step 3: obtain a long-range potential
  BreakBase* breaker = create_break(box, doc);

  // step 4: construct finite size correction integrals
  int nrule = 4;
  int nk = 64;
  RealType kmax = breaker->get_kc();
  RealType dk = kmax/nk;
  vector<RealType> kmags(nk);
  for (int ik=0; ik<nk; ik+=1)
  {
    kmags[ik] = ik*dk;
  }
  //vector<RealType> intvals = spherical_integral(boxspl3d, kmags, nrule);
  ofstream ofs;
  //ofs.open("avesk.dat", ofstream::out);
  ofs.open("vk.dat");
  for (int ik=0; ik<nk; ik+=1)
  {
    //ofs << kmags[ik] << " " << intvals[ik] << endl;
    ofs << kmags[ik] << " " << breaker->evaluate_fklr(kmags[ik]) << endl;
  }
  ofs.close();

  OHMMS::Controller->finalize();
}
