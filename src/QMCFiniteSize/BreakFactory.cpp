#include "QMCFiniteSize/BreakFactory.h"
namespace qmcplusplus
{
BreakBase* BreakFactory::create_break(
  Uniform3DGridLayout box,
  xmlXPathContextPtr doc
)
{
  string basis = "esler";
  int nknot;
  RealType rc, kc, kcut, kmax;
  // set defaults
  box.SetLRCutoffs();
  rc = box.LR_rc;
  kc = box.LR_kc;
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
} // qmcplusplus
