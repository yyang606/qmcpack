#ifndef QMCPLUSPLUS_ESLER_BREAK_H
#define QMCPLUSPLUS_ESLER_BREAK_H
#include "LongRange/LPQHIBasis.h"
#include "QMCFiniteSize/BreakBase.h"
namespace qmcplusplus
{

typedef LPQHIBasis                                        BasisType;
typedef std::function< RealType(RealType, RealType) > EslerFuncType;

class EslerBreak
{
 public:
  EslerBreak(
    EslerFuncType fxk,
    BoxType box,
    BreakSpec params
  );
  ~EslerBreak();
  RealType get_chisq(){return chisq_;};
  RealType get_rc(){return params_.rc;};
  RealType get_kc(){return params_.kc;};

  // goal in life: evaluate long-range potential at k
  RealType evaluate_fklr(RealType k);
 private:
  EslerFuncType fxk_;
  BasisType basis_;
  LRBreakup<BasisType>* handler_;
  BoxType box_;  // define reciprocal lattice
  std::vector<RealType> coefs_;
  RealType chisq_;
  BreakSpec params_;
  friend std::ostream& operator<<(std::ostream& os, const EslerBreak& breaker);
};

}
#endif
