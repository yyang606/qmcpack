#include "QMCFiniteSize/lr_routines.h"
#include "LongRange/KContainer.h"

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
  app_log() << endl;
  breaker->report(app_log());
  app_log() << "  chi^2  = " << scientific << breaker->get_chisq() << endl;
  app_log() << setprecision(10);

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
  vector<RealType> intvals = spherical_integral(boxspl3d, kmags, nrule);
  ofstream ofs, ofv, ofi;
  ofs.open("avesk.dat");
  ofv.open("vk.dat");
  ofi.open("vint.dat");
  for (int ik=0; ik<nk; ik+=1)
  {
    RealType sk = intvals[ik];
    RealType vklr = breaker->evaluate_fklr(kmags[ik]);
    RealType norm = box.Volume/(2*M_PI*M_PI);
    ofs << kmags[ik] << " " << sk << endl;
    ofv << kmags[ik] << " " << vklr << endl;
    ofi << kmags[ik] << " " << norm*0.5*std::pow(kmags[ik], 2)*vklr*sk << endl;
  }
  ofs.close();
  ofv.close();
  ofi.close();

  // step 5: dump finite size sum
  KContainer kvecs;
  kvecs.UpdateKLists(box, breaker->get_kc());
  nk = kvecs.ksq.size();
  kmags.resize(nk);
  for (int ik=0; ik<nk; ik++)
  {
    kmags[ik] = sqrt(kvecs.ksq[ik]);
  }
  intvals = spherical_integral(boxspl3d, kmags, nrule);
  RealType vsum = 0.0;
  for (int ik=0; ik<nk; ik++)
  {
    RealType kmag, vklr, sk, uk, uklr;
    kmag = kmags[ik];
    vklr = breaker->evaluate_fklr(kmag);
    sk = intvals[ik];
    vsum += 0.5*vklr*sk;
  }
  vsum = vsum;
  app_log() << "  vsum = " << vsum << endl;

  OHMMS::Controller->finalize();
}
