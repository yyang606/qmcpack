#include "QMCFiniteSize/BreakFactory.h"
namespace qmcplusplus
{
BreakBase* BreakFactory::create_break(
  Uniform3DGridLayout box,
  xmlXPathContextPtr doc
)
{
  string basis;
  int nknot;
  RealType rc, kc, kcut, kmax;

  // set defaults
  basis = "esler";
  box.SetLRCutoffs();
  rc = box.LR_rc;
  //kc = box.LR_kc; // !!!! dangerous default if S(k) is on rectangular lattice
  // spline3d needs to extrapolate if max|k| is larger along x than along y
  RealType kc0 = default_kc(box);
  kc = kc0;
  nknot = 15;
  kcut = 60*M_PI*std::pow((box.Volume), -1./3);
  kmax = 6000./rc;

  // parse user input
  xmlNodePtr node;
  node = find("//lrbreak/short-range/basis", doc);
  if (node) putContent(basis, node);
  node = find("//lrbreak/short-range/rc", doc);
  if (node) putContent(rc, node);
  node = find("//lrbreak/short-range/nknot", doc);
  if (node) putContent(nknot, node);
  node = find("//lrbreak/long-range/kc", doc);
  if (node) putContent(kc, node);
  node = find("//lrbreak/long-range/kcut", doc);
  if (node) putContent(kcut, node);
  node = find("//lrbreak/long-range/kmax", doc);
  if (node) putContent(kmax, node);

  // warn user about potentially bad input (comment out to continue)
  if (kc > kc0)
  {
    APP_ABORT("kc = " << kc << " > default "  <<  kc0 << endl << " Note: "
      << "box kc may require 3D S(k) spline to extrapolate in non-cubic cell"
      << endl);
   }
  if (kcut < kc) APP_ABORT("continuum kcut = " << kcut << " < kc = " << kc);
  if (kmax < kcut) APP_ABORT("kmax < kcut");

  // instantiate potential breaker
  BreakSpec params;
  params.basis = basis;
  params.rc = rc;
  params.kc = kc;
  params.kcut = kcut;
  params.kmax = kmax;
  params.nknot = nknot;
  BreakBase* breaker(0);
  if (basis == "natoli")
  {
    NatoliCoul fk(box.Volume);
    breaker = new NatoliBreak(fk, box, params);
  } else
  if (basis == "esler")
  {
    EslerCoul fxk(box.Volume);
    breaker = new EslerBreak(fxk, box, params);
  } else
  {
    APP_ABORT("unknown basis " << basis << endl);
  }
  return breaker;
}
RealType BreakFactory::default_kc(Uniform3DGridLayout box)
{ // reduce LR_kc of vklr so that S(k) 3D spline never has to extrapolate
  RealType kc=box.LR_kc;
  RealType max_blat = 0.0;
  for (int idim=0; idim<3; idim++)
  {
    RealType blat = 2*M_PI*std::sqrt( dot(box.Gv[idim], box.Gv[idim]) );
    if (blat>max_blat) max_blat = blat;
  }
  return kc-max_blat;
}
} // qmcplusplus
