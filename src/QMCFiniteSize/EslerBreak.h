#ifndef QMCPLUSPLUS_ESLER_BREAK_H
#define QMCPLUSPLUS_ESLER_BREAK_H
#include "LongRange/LPQHIBasis.h"
#include "QMCFiniteSize/BreakBase.h"
namespace qmcplusplus
{

typedef LPQHIBasis                                        BasisType;
typedef std::function< RealType(RealType, RealType) > EslerFuncType;

class EslerBreak : public BreakBase
{
 public:
  EslerBreak(
    EslerFuncType fxk,
    BoxType box,
    BreakSpec params
  );
  ~EslerBreak();
  void report(std::ostream& os);

  // goal in life: evaluate long-range potential at k
  RealType evaluate_fklr(RealType k);
 private:
  EslerFuncType fxk_;
  BasisType basis_;
  LRBreakup<BasisType>* handler_;
  BoxType box_;  // define reciprocal lattice
  std::vector<RealType> coefs_;
};

}
#endif
