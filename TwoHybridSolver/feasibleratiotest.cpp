//
// File: feasibleratiotest.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

// Include Files
#include "feasibleratiotest.h"
#include "rt_nonfinite.h"
#include "xnrm2.h"
#include <cmath>
#include <string.h>

// Function Definitions
//
// Arguments    : const double solution_xstar[5]
//                const double solution_searchDir[5]
//                int workingset_nVar
//                const double workingset_lb[5]
//                const double workingset_ub[5]
//                const int workingset_indexLB[5]
//                const int workingset_indexUB[5]
//                const int workingset_sizes[5]
//                const int workingset_isActiveIdx[6]
//                const boolean_T workingset_isActiveConstr[9]
//                const int workingset_nWConstr[5]
//                boolean_T isPhaseOne
//                double *alpha
//                boolean_T *newBlocking
//                int *constrType
//                int *constrIdx
// Return Type  : void
//
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
    boolean_T *newBlocking, int *constrType, int *constrIdx)
{
  double denomTol;
  double phaseOneCorrectionP;
  double phaseOneCorrectionX;
  double pk_corrected;
  double ratio;
  *alpha = 1.0E+30;
  *newBlocking = false;
  *constrType = 0;
  *constrIdx = 0;
  denomTol = 2.2204460492503131E-13 *
             internal::blas::xnrm2(workingset_nVar, solution_searchDir);
  if (workingset_nWConstr[3] < workingset_sizes[3]) {
    int i;
    phaseOneCorrectionX =
        static_cast<double>(isPhaseOne) * solution_xstar[workingset_nVar - 1];
    phaseOneCorrectionP = static_cast<double>(isPhaseOne) *
                          solution_searchDir[workingset_nVar - 1];
    i = workingset_sizes[3];
    for (int idx{0}; idx <= i - 2; idx++) {
      int i1;
      i1 = workingset_indexLB[idx];
      pk_corrected = -solution_searchDir[i1 - 1] - phaseOneCorrectionP;
      if ((pk_corrected > denomTol) &&
          (!workingset_isActiveConstr[(workingset_isActiveIdx[3] + idx) - 1])) {
        ratio = (-solution_xstar[i1 - 1] - workingset_lb[i1 - 1]) -
                phaseOneCorrectionX;
        pk_corrected =
            std::fmin(std::abs(ratio), 1.0E-6 - ratio) / pk_corrected;
        if (pk_corrected < *alpha) {
          *alpha = pk_corrected;
          *constrType = 4;
          *constrIdx = idx + 1;
          *newBlocking = true;
        }
      }
    }
    i = workingset_indexLB[workingset_sizes[3] - 1] - 1;
    pk_corrected = -solution_searchDir[i];
    if ((pk_corrected > denomTol) &&
        (!workingset_isActiveConstr
             [(workingset_isActiveIdx[3] + workingset_sizes[3]) - 2])) {
      ratio = -solution_xstar[i] - workingset_lb[i];
      pk_corrected = std::fmin(std::abs(ratio), 1.0E-6 - ratio) / pk_corrected;
      if (pk_corrected < *alpha) {
        *alpha = pk_corrected;
        *constrType = 4;
        *constrIdx = workingset_sizes[3];
        *newBlocking = true;
      }
    }
  }
  if (workingset_nWConstr[4] < 2) {
    phaseOneCorrectionX =
        static_cast<double>(isPhaseOne) * solution_xstar[workingset_nVar - 1];
    phaseOneCorrectionP = static_cast<double>(isPhaseOne) *
                          solution_searchDir[workingset_nVar - 1];
    pk_corrected =
        solution_searchDir[workingset_indexUB[0] - 1] - phaseOneCorrectionP;
    if ((pk_corrected > denomTol) &&
        (!workingset_isActiveConstr[workingset_isActiveIdx[4] - 1])) {
      ratio = (solution_xstar[workingset_indexUB[0] - 1] -
               workingset_ub[workingset_indexUB[0] - 1]) -
              phaseOneCorrectionX;
      pk_corrected = std::fmin(std::abs(ratio), 1.0E-6 - ratio) / pk_corrected;
      if (pk_corrected < *alpha) {
        *alpha = pk_corrected;
        *constrType = 5;
        *constrIdx = 1;
        *newBlocking = true;
      }
    }
    pk_corrected =
        solution_searchDir[workingset_indexUB[1] - 1] - phaseOneCorrectionP;
    if ((pk_corrected > denomTol) &&
        (!workingset_isActiveConstr[workingset_isActiveIdx[4]])) {
      ratio = (solution_xstar[workingset_indexUB[1] - 1] -
               workingset_ub[workingset_indexUB[1] - 1]) -
              phaseOneCorrectionX;
      pk_corrected = std::fmin(std::abs(ratio), 1.0E-6 - ratio) / pk_corrected;
      if (pk_corrected < *alpha) {
        *alpha = pk_corrected;
        *constrType = 5;
        *constrIdx = 2;
        *newBlocking = true;
      }
    }
  }
  if (!isPhaseOne) {
    if ((*newBlocking) && (*alpha > 1.0)) {
      *newBlocking = false;
    }
    *alpha = std::fmin(*alpha, 1.0);
  }
}

} // namespace qpactiveset
} // namespace coder
} // namespace optim
} // namespace coder

//
// File trailer for feasibleratiotest.cpp
//
// [EOF]
//
