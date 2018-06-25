#include "QMCFiniteSize/fsc_routines.h"
#include "QMCFiniteSize/Spline3DFactory.h"
#include "QMCFiniteSize/SphericalAverage3D.h"
#include "QMCFiniteSize/BreakFactory.h"
#include "QMCFiniteSize/Quad1D.h"
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

  // step 3: obtain a long-range potential
  BreakFactory break_factory;
  BreakBase* breaker = break_factory.create_break(box, doc);
  app_log() << endl;
  breaker->report(app_log());
  app_log() << "  chi^2  = " << scientific << breaker->get_chisq() << endl;
  app_log() << setprecision(10);

  // step 4: setup integral parameters
  int nrule=4, nk=64;
  xmlNodePtr node;
  node = find("//output_grid/nk", doc);
  if (node) putContent(nk, node);
  node = find("//spherical_average/nrule", doc);
  if (node) putContent(nrule, node);
  app_log() << " nrule = " << nrule << endl;
  SphericalAverage3D sphavg(nrule);
  RealType kmax = box.LR_kc;

  // step 5: do finite size correction integrals
  RealType vint = 0;
  Quad1D quad1d(0, kmax, nk);
  // output useful stuff
  ofstream ofs, ofv, ofi;
  ofs.open("avesk.dat");
  ofv.open("vk.dat");
  ofi.open("vint.dat");
  for (int ik=0; ik<nk; ik+=1)
  {
    RealType kmag = quad1d.x[ik];
    RealType sk = sphavg(*boxspl3d, kmag);
    RealType vklr = breaker->evaluate_fklr(kmag);
    ofs << kmag << " " << sk << endl;
    ofv << kmag << " " << vklr << endl;
    RealType val = 0.5*sk*breaker->evaluate_fklr(kmag);
    ofi << kmag << " " << val << endl;
    //vint += 0.5*vklr*sk*quad1d.w[ik];
    vint += val*quad1d.w[ik];
  }
  ofs.close();
  ofv.close();
  ofi.close();
  app_log() << " vint = " << vint << endl;

  OHMMS::Controller->finalize();
}
