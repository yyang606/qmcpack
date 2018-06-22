#ifndef QMCPLUSPLUS_ESLER_BREAK_H
#define QMCPLUSPLUS_ESLER_BREAK_H
#include <vector>
#include <functional>
#include <ostream>
#include "LongRange/LPQHIBasis.h"
#include "LongRange/LRBreakup.h"
namespace qmcplusplus
{

typedef LPQHIBasis                                   BasisType;
typedef Uniform3DGridLayout                            BoxType;
typedef QMCTraits::RealType                           RealType;
typedef std::function< RealType(RealType, RealType) > FuncType;

struct BreakSpec
{
  int nknot;
  double rc, kc, kcut, kmax;
};

class EslerBreak
{
 public:
  EslerBreak(
    FuncType fxk,
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
  BasisType basis_;
  BoxType box_;  // define reciprocal lattice
  LRBreakup<BasisType>* handler_;
  std::vector<RealType> coefs_;
  RealType chisq_;
  FuncType fxk_;
  friend std::ostream& operator<<(std::ostream& os, const EslerBreak& breaker);
  BreakSpec params_;
};

}
#endif
