#include "QMCFiniteSize/EslerBreak.h"

using namespace std;

namespace qmcplusplus
{
// fxk(k, rc)
EslerBreak::EslerBreak(
  FuncType fxk,
  BoxType box,
  BreakSpec params
) : fxk_(fxk), basis_(box), box_(box), params_(params), chisq_(-1)
{
  // define short-range basis, coefs_ still need to be filled
  basis_.set_Lattice(box);
  basis_.set_NumKnots(params.nknot);
  basis_.set_rc(params.rc);
  coefs_.resize(basis_.NumBasisElem());
  // define long-range basis i.e. kvectors
  handler_ = new LRBreakup<BasisType>(basis_);
  handler_->SetupKVecs(params.kc, params.kcut, params.kmax);
  // solve for short-range coefficients
  int nk = handler_->KList.size();
  vector<RealType> xk(nk);
  for (int ik=0; ik<nk; ik++)
  { // KList has 2 columns: kshell magnitude, kshell weight
    xk[ik] = fxk(handler_->KList[ik][0], basis_.get_rc());
  }
  chisq_ = handler_->DoBreakup(xk.data(), coefs_.data());
}

EslerBreak::~EslerBreak()
{
  delete handler_;
}

RealType EslerBreak::evaluate_fklr(RealType k)
{
  RealType val = fxk_(k, basis_.get_rc());
  for (int n=0; n<basis_.NumBasisElem(); n++)
  {
    val -= coefs_[n]*basis_.c(n, k);
  }
  return val;
}

ostream& operator<<(ostream& os, const EslerBreak& breaker)
{
  os << " Esler long-range breakup" << endl;
  os << " ------------------------" << endl;
  os << "  nknot  = " << breaker.params_.nknot << endl;
  os << "  rc     = " << breaker.params_.rc    << endl;
  os << "  kc     = " << breaker.params_.kc    << endl;
  os << "  kcut   = " << breaker.params_.kcut  << endl;
  os << "  kmax   = " << breaker.params_.kmax  << endl;
  return os;
}
}
