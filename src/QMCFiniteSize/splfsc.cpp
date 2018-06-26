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
  RealType kmax = breaker->get_kc();

  /*
  // !!!! HACK: use bare Coulomb below kmin
  RealType kmin = 0.05;
  app_log() << " !!!! HACK: use bare Coulomb for k<kmin" << endl;

  // step 5: do finite size correction integrals
  RealType norm = box.Volume/(2*M_PI*M_PI);
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

  kmin = -kmax;
  ofstream ofx, ofy, ofz;
  ofx.open("sx.dat");
  ofy.open("sy.dat");
  ofz.open("sz.dat");
  RealType dk = (kmax-kmin)/(nk-1);
  for (int ik=0; ik<nk; ik++)
  {
    RealType kmag = kmin+ik*dk;
    ofx << kmag << " " << (*boxspl3d)(kmag, 0, 0) << endl;
    ofy << kmag << " " << (*boxspl3d)(0, kmag, 0) << endl;
    ofz << kmag << " " << (*boxspl3d)(0, 0, kmag) << endl;
  }
  ofx.close();
  ofy.close();
  ofz.close();
  */

  // step 5: do finite size correction sums
  KContainer kvecs;
  kvecs.UpdateKLists(box, kmax);
  RealType vsum = 0.0;
  for (int ik=0; ik<kvecs.ksq.size(); ik++)
  {
    RealType kmag = std::sqrt(kvecs.ksq[ik]);
    RealType vklr = breaker->evaluate_fklr(kmag);
    RealType sk = sphavg(*boxspl3d, kmag);
    vsum += 0.5*vklr*sk;
  }

  // step 6: evaluate spline on kshells, check sums
  ofstream ofg;
  ofg.open("integrand.dat");
  int nks = kvecs.kshell.size()-1;  // number of kshells
  RealType vsum1 = 0.0;
  // +1 for left value boundary condition
  vector<RealType> kmags(nks+1, 0);
  vector<RealType> integrand1d(nks+1, 0);
  for (int iks=0; iks<nks; iks++)
  {
    RealType kmag = std::sqrt(kvecs.ksq[kvecs.kshell[iks]]);
    RealType vklr = breaker->evaluate_fklr(kmag);
    RealType sk = sphavg(*boxspl3d, kmag);
    RealType val = 0.5*vklr*sk;
    for (int ik=kvecs.kshell[iks]; ik<kvecs.kshell[iks+1]; ik++)
    {
      vsum1 += val;
    }
    kmags[iks+1] = kmag;
    integrand1d[iks+1] = std::pow(kmag, 2)*val;
    ofg << kmags[iks] << " " << integrand1d[iks] << endl;
  }
  ofg.close();
  if (abs(vsum1-vsum)>1e-16) APP_ABORT("vsum mismatch");

  Spline1D spline1d(kmags, integrand1d);

  nk = 128;

  RealType norm = box.Volume/(2*M_PI*M_PI);
  Quad1D quad1d(0, kmax, nk);
  RealType vint = 0.0;
  for (int ik=0; ik<nk; ik+=1)
  {
    RealType kmag = quad1d.x[ik];
    vint += norm*spline1d(kmag)*quad1d.w[ik];
  }

  app_log() << " vint = " << vint << endl;
  app_log() << " vsum  = " << vsum << endl;
  app_log() << " dvlr = " << vint - vsum << endl;

  OHMMS::Controller->finalize();
}
