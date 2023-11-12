//
// File: checkStoppingAndUpdateFval.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

// Include Files
#include "checkStoppingAndUpdateFval.h"
#include "TwoHybridSolver_internal_types.h"
#include "computeFval_ReuseHx.h"
#include "feasibleX0ForWorkingSet.h"
#include "rt_nonfinite.h"
#include <algorithm>
#include <cmath>
#include <string.h>

// Function Definitions
//
// Arguments    : int *activeSetChangeID
//                const double f[5]
//                i_struct_T *solution
//                d_struct_T *memspace
//                const struct_T *objective
//                const j_struct_T *workingset
//                e_struct_T *qrmanager
//                double options_ObjectiveLimit
//                int runTimeOptions_MaxIterations
//                boolean_T updateFval
// Return Type  : void
//
namespace coder {
namespace optim {
namespace coder {
namespace qpactiveset {
namespace stopping {
void checkStoppingAndUpdateFval(int *activeSetChangeID, const double f[5],
                                i_struct_T *solution, d_struct_T *memspace,
                                const struct_T *objective,
                                const j_struct_T *workingset,
                                e_struct_T *qrmanager,
                                double options_ObjectiveLimit,
                                int runTimeOptions_MaxIterations,
                                boolean_T updateFval)
{
  int nVar;
  solution->iterations++;
  nVar = objective->nvar - 1;
  if ((solution->iterations >= runTimeOptions_MaxIterations) &&
      ((solution->state != 1) || (objective->objtype == 5))) {
    solution->state = 0;
  }
  if (solution->iterations - solution->iterations / 50 * 50 == 0) {
    double x;
    int b_idxUB_tmp;
    int idx;
    int idxLB;
    int idxUB_tmp;
    int mLB;
    mLB = workingset->sizes[3];
    x = 0.0;
    for (idx = 0; idx < mLB; idx++) {
      idxLB = workingset->indexLB[idx] - 1;
      x = std::fmax(x, -solution->xstar[idxLB] - workingset->lb[idxLB]);
    }
    idxUB_tmp = workingset->indexUB[0] - 1;
    x = std::fmax(x, solution->xstar[idxUB_tmp] - workingset->ub[idxUB_tmp]);
    b_idxUB_tmp = workingset->indexUB[1] - 1;
    x = std::fmax(x,
                  solution->xstar[b_idxUB_tmp] - workingset->ub[b_idxUB_tmp]);
    solution->maxConstr = x;
    if (x > 1.0E-6) {
      double constrViolation_new;
      boolean_T nonDegenerateWset;
      if (0 <= nVar) {
        std::copy(&solution->xstar[0], &solution->xstar[nVar + 1],
                  &solution->searchDir[0]);
      }
      nonDegenerateWset = initialize::feasibleX0ForWorkingSet(
          memspace->workspace_double, solution->searchDir, workingset,
          qrmanager);
      if ((!nonDegenerateWset) && (solution->state != 0)) {
        solution->state = -2;
      }
      *activeSetChangeID = 0;
      mLB = workingset->sizes[3];
      constrViolation_new = 0.0;
      for (idx = 0; idx < mLB; idx++) {
        idxLB = workingset->indexLB[idx] - 1;
        constrViolation_new =
            std::fmax(constrViolation_new,
                      -solution->searchDir[idxLB] - workingset->lb[idxLB]);
      }
      constrViolation_new =
          std::fmax(constrViolation_new,
                    solution->searchDir[idxUB_tmp] - workingset->ub[idxUB_tmp]);
      constrViolation_new =
          std::fmax(constrViolation_new, solution->searchDir[b_idxUB_tmp] -
                                             workingset->ub[b_idxUB_tmp]);
      if (constrViolation_new < x) {
        if (0 <= nVar) {
          std::copy(&solution->searchDir[0], &solution->searchDir[nVar + 1],
                    &solution->xstar[0]);
        }
        solution->maxConstr = constrViolation_new;
      }
    }
  }
  if ((options_ObjectiveLimit > rtMinusInf) && updateFval) {
    solution->fstar = Objective::computeFval_ReuseHx(
        objective, memspace->workspace_double, f, solution->xstar);
    if ((solution->fstar < options_ObjectiveLimit) &&
        ((solution->state != 0) || (objective->objtype != 5))) {
      solution->state = 2;
    }
  }
}

} // namespace stopping
} // namespace qpactiveset
} // namespace coder
} // namespace optim
} // namespace coder

//
// File trailer for checkStoppingAndUpdateFval.cpp
//
// [EOF]
//
