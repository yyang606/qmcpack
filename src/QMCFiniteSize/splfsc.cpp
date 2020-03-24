#include "QMCFiniteSize/fsc_routines.h"
#include "QMCFiniteSize/Spline3DFactory.h"
#include "QMCFiniteSize/SphericalAverage3D.h"
#include "QMCFiniteSize/BreakFactory.h"
#include "QMCFiniteSize/Quad1D.h"
#include "QMCFiniteSize/Spline1D.h"
#include "einspline/nubspline_create.h"
#include "einspline/nubspline_eval_d.h"
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
  xmlNodePtr sc_node = find("//simulationcell", doc);
  Uniform3DGridLayout box = create_box(sc_node);

  // step 2: obtain spline on regular grid (in reciprocal lattice units)
  vector<vector<RealType>> mat = loadtxt(fdat_name);
  Spline3DFactory spline_factory;
  Ugrid3D grid3d = spline_factory.create_ugrid3d(doc, "input_grid");
  NaturalSpline3DInBox* boxspl3d = spline_factory.create_boxspl3d(
    grid3d, mat, box);

  // step 3: obtain a long-range potential
  BreakFactory break_factory;
  BreakBase* breaker = break_factory.create_break(box, doc);
  breaker->report(app_log());
  RealType kmax = breaker->get_kc();

  // step 4: setup integral parameters
  int nrule=4, nk=64;
  xmlNodePtr node;
  node = find("//output_grid/nk", doc);
  if (node) putContent(nk, node);
  node = find("//spherical_average/nrule", doc);
  if (node) putContent(nrule, node);
  app_log() << " nrule = " << nrule << endl;
  SphericalAverage3D sphavg(nrule);
  RealType dvlr;

  // step 5: do finite size sum, store spherical average at kshells
  //  the spline is most accurate at kshells
  KContainer kvecs;
  kvecs.UpdateKLists(box, kmax);
  RealType vsum = 0.0;
  int nks = kvecs.kshell.size()-1;  // number of kshells
  // +1 for left value boundary condition S(k->0) = 0
  vector<RealType> kmags(nks+1, 0);
  vector<RealType> vint1d(nks+1, 0);
  for (int iks=0; iks<nks; iks++)
  {
    RealType kmag = std::sqrt(kvecs.ksq[kvecs.kshell[iks]]);
    RealType vklr = breaker->evaluate_fklr(kmag);
    RealType sk = sphavg(*boxspl3d, kmag);
    RealType val = 0.5*vklr*sk;
    for (int ik=kvecs.kshell[iks]; ik<kvecs.kshell[iks+1]; ik++)
    { // KContainer should have a list of pre-calculated weights
      vsum += val;
    }
    kmags[iks+1] = kmag;
    vint1d[iks+1] = std::pow(kmag, 2)*val;
  }

  // step 6: estimate thermodynamic limit of the sum
  //  !!!! this is the most tricky step. Try a few approaches
  RealType vint = 0.0;
  RealType norm = box.Volume/(2*M_PI*M_PI);

  // attempt 1: directly use 3D spline
  // problem: break->evaluate_vklr is numerically unstable for k < 0.05
  //
  // !!!! HACK: use bare Coulomb below kmin
  RealType kmin = 0.05;
  app_log() << " !!!! HACK: use bare Coulomb for k<kmin" << endl;

  Quad1D quad1d(0, kmax, nk);
  vint = 0.0;
  ofstream ofs, ofv, ofi;
  ofs.open("avesk.dat");
  ofv.open("vklr.dat");
  ofi.open("vint.dat");
  for (int ik=0; ik<nk; ik+=1)
  {
    RealType kmag = quad1d.x[ik];
    RealType sk = sphavg(*boxspl3d, kmag);
    RealType vklr = breaker->evaluate_fklr(kmag);
    if (kmag < kmin) vklr = 4*M_PI/(kmag*kmag)/box.Volume;
    ofs << kmag << " " << sk << endl;
    ofv << kmag << " " << vklr << endl;
    RealType val = 0.5*std::pow(kmag, 2)*sk*vklr*norm;
    ofi << kmag << " " << val << endl;
    vint += val*quad1d.w[ik];
  }
  ofs.close();
  ofv.close();
  ofi.close();
  dvlr = vint - vsum;
  app_log() << endl;
  app_log() << " Directly using 3D spline " << endl;
  app_log() << " ========================" << endl;
  app_log() << " vint = " << vint << endl;
  app_log() << " vsum = " << vsum << endl;
  app_log() << " dvlr = " << dvlr << endl;
  app_log() << endl;

  // attempt 2: respline integrand in 1D
  // problem: spline quality depends on the shape of the integrand as k->0
  //  monotonic = good, oscilitory = bad!
  //
  ofstream ofg;
  ofg.open("vint1d.dat");
  Spline1D spline1d(kmags, vint1d);
  for (int iks=0; iks<vint1d.size(); iks++)
  {
    ofg << kmags[iks] << " " << vint1d[iks]*norm << endl;
  }
  ofg.close();

  //Quad1D quad1d(0, kmags[nks], nk);
  vint = 0.0;
  ofg.open("fvint1d.dat");
  for (int ik=0; ik<nk; ik+=1)
  {
    RealType kmag = quad1d.x[ik];
    RealType val = spline1d(kmag)*norm;
    ofg << kmag << " " << val << endl;
    vint += val*quad1d.w[ik];
  }
  ofg.close();

  dvlr = vint - vsum;
  app_log() << endl;
  app_log() << " Respline integrand in 1D " << endl;
  app_log() << " =========================" << endl;
  app_log() << " vint = " << vint << endl;
  app_log() << " vsum = " << vsum << endl;
  app_log() << " dvlr = " << dvlr << endl;
  app_log() << endl;

  // attempt 3: do sum in a bigger box
  node = find("//bigcell", doc);
  Uniform3DGridLayout bigbox = create_box(node);
  KContainer bigkvecs;
  bigkvecs.UpdateKLists(bigbox, kmax);
  vint = 0.0;
  for (int ik=0; ik<bigkvecs.ksq.size(); ik++)
  {
    RealType kmag = std::sqrt(bigkvecs.ksq[ik]);
    RealType vklr = breaker->evaluate_fklr(kmag)*box.Volume/bigbox.Volume;
    RealType sk = sphavg(*boxspl3d, kmag);
    vint += 0.5*vklr*sk;
  }

  app_log() << endl;
  app_log() << " Do bigger cell, then 1/N extrap." << endl;
  app_log() << " ================================" << endl;
  int nvol = bigbox.Volume/box.Volume;
  dvlr = (nvol*vint-vsum)/(nvol-1.) - vsum;
  app_log() << " vbgs = " << vint << endl;
  app_log() << " vsum = " << vsum << endl;
  app_log() << " dvlr = " << dvlr << endl;
  app_log() << endl;

  OHMMS::Controller->finalize();
}
