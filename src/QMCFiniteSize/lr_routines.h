#include "QMCFiniteSize/KspaceFunctions.hpp"
#include "QMCFiniteSize/EslerBreak.h"
#include "QMCFiniteSize/NatoliBreak.h"
#include "QMCFiniteSize/fsc_routines.h"
#ifndef LR_COMMON_H
#define LR_COMMON_H
using namespace qmcplusplus;
NatoliBreak create_natoli_break(
  Uniform3DGridLayout box,
  xmlXPathContextPtr doc
);
#endif
