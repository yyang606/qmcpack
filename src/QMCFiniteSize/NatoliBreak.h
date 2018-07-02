#ifndef QMCPLUSPLUS_NATOLI_BREAK_H
#define QMCPLUSPLUS_NATOLI_BREAK_H
#include "LongRange/LPQHISRCoulombBasis.h"
#include "QMCFiniteSize/BreakBase.h"
namespace qmcplusplus
{
typedef LPQHISRCoulombBasis                 NatoliBasisType;
typedef std::function< RealType(RealType) >        FuncType;

class NatoliBreak : public BreakBase
{
 public:
  NatoliBreak(
    FuncType fk,
    BoxType box,
    BreakSpec params
  );
  ~NatoliBreak();
  void report(std::ostream& os);

  // goal in life: evaluate long-range potential at k
  RealType evaluate_fklr(RealType k);

 private:
  FuncType fk_;
  NatoliBasisType basis_;
  LRBreakup<NatoliBasisType>* handler_;
  BoxType box_;
  std::vector<RealType> coefs_;
};
} // qmcplusplus
#endif