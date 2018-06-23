#include "QMCFiniteSize/EslerBreak.h"

using namespace std;

namespace qmcplusplus
{
EslerBreak::EslerBreak(
  EslerFuncType fxk,
  BoxType box,
  BreakSpec params
) : fxk_(fxk), basis_(box), box_(box), BreakBase(params)
{
  // define short-range basis, coefs_ still needs to be filled
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
  RealType val = fxk_(k, basis_.get_rc())-basis_.fk(k, coefs_);
  return val;
}

void EslerBreak::report(ostream& os)
{
  os << " Esler long-range breakup" << endl;
  os << " ------------------------" << endl;
  os << params_ << endl;
}
}
