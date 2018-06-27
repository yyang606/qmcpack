#include "QMCFiniteSize/NatoliBreak.h"

using namespace std;

namespace qmcplusplus
{
NatoliBreak::NatoliBreak(
  FuncType fk,
  BoxType box,
  BreakSpec params
) : fk_(fk), basis_(box), box_(box), BreakBase(params)
{
  // define short-range basis, coefs_ still needs to be filled
  basis_.set_Lattice(box);
  basis_.set_NumKnots(params.nknot);
  basis_.set_rc(params.rc);
  coefs_.resize(basis_.NumBasisElem());
  // define long-range basis i.e. kvectors
  handler_ = new LRBreakup<NatoliBasisType>(basis_);
  handler_->SetupKVecs(params.kc, params.kcut, params.kmax);
  // solve for short-range coefficients
  int nk = handler_->KList.size();
  vector<RealType> vk(nk);
  for (int ik=0; ik<nk; ik++)
  { // KList has 2 columns: kshell magnitude, kshell weight
    vk[ik] = fk(handler_->KList[ik][0]);
  }
  // fit constraints
  int nbasis = basis_.NumBasisElem();
  vector<RealType> constraints(nbasis, 1.0);
  constraints[0] = 0;
  constraints[1] = 0;
  constraints[2] = 0;
  constraints[nbasis-1] = 0;
  constraints[nbasis-2] = 0;
  constraints[nbasis-3] = 0;
  coefs_[0] = 1.0;
  coefs_[1] = 0.0;
  coefs_[2] = 0.0;
  coefs_[nbasis-1] = 0.0;
  coefs_[nbasis-2] = 0.0;
  coefs_[nbasis-3] = 0.0;
  chisq_ = handler_->DoBreakup(vk.data(), coefs_.data(), constraints.data());
}
NatoliBreak::~NatoliBreak()
{
  delete handler_;
}
RealType NatoliBreak::evaluate_fklr(RealType k)
{
  RealType val = fk_(k)-basis_.fk(k, coefs_);
  return val;
}
} // qmcplusplus
