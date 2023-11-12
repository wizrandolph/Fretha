//
// File: feasibleratiotest.h
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

#ifndef FEASIBLERATIOTEST_H
#define FEASIBLERATIOTEST_H

// Include Files
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
namespace optim {
namespace coder {
namespace qpactiveset {
void feasibleratiotest(
    const double solution_xstar[5], const double solution_searchDir[5],
    int workingset_nVar, const double workingset_lb[5],
    const double workingset_ub[5], const int workingset_indexLB[5],
    const int workingset_indexUB[5], const int workingset_sizes[5],
    const int workingset_isActiveIdx[6],
    const boolean_T workingset_isActiveConstr[9],
    const int workingset_nWConstr[5], boolean_T isPhaseOne, double *alpha,
    boolean_T *newBlocking, int *constrType, int *constrIdx);

}
} // namespace coder
} // namespace optim
} // namespace coder

#endif
//
// File trailer for feasibleratiotest.h
//
// [EOF]
//
