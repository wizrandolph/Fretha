//
// File: driver1.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

// Include Files
#include "driver1.h"
#include "PresolveWorkingSet.h"
#include "TwoHybridSolver_internal_types.h"
#include "computeFval.h"
#include "iterate.h"
#include "rt_nonfinite.h"
#include "setProblemType.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <string.h>

// Function Definitions
//
// Arguments    : const double H[16]
//                const double f[5]
//                i_struct_T *solution
//                d_struct_T *memspace
//                j_struct_T *workingset
//                e_struct_T *qrmanager
//                f_struct_T *cholmanager
//                struct_T *objective
//                g_struct_T *options
//                int runTimeOptions_MaxIterations
// Return Type  : void
//
namespace coder {
namespace optim {
namespace coder {
namespace qpactiveset {
void driver(const double H[16], const double f[5], i_struct_T *solution,
            d_struct_T *memspace, j_struct_T *workingset, e_struct_T *qrmanager,
            f_struct_T *cholmanager, struct_T *objective, g_struct_T *options,
            int runTimeOptions_MaxIterations)
{
  int idx;
  int idxLB;
  int nVar;
  boolean_T guard1{false};
  solution->iterations = 0;
  nVar = workingset->nVar - 1;
  guard1 = false;
  if (workingset->probType == 3) {
    idxLB = workingset->sizes[3];
    for (idx = 0; idx < idxLB; idx++) {
      if (workingset->isActiveConstr[idx]) {
        solution->xstar[workingset->indexLB[idx] - 1] =
            -workingset->lb[workingset->indexLB[idx] - 1];
      }
    }
    if (workingset->isActiveConstr[workingset->isActiveIdx[4] - 1]) {
      solution->xstar[workingset->indexUB[0] - 1] =
          workingset->ub[workingset->indexUB[0] - 1];
    }
    if (workingset->isActiveConstr[workingset->isActiveIdx[4]]) {
      solution->xstar[workingset->indexUB[1] - 1] =
          workingset->ub[workingset->indexUB[1] - 1];
    }
    initialize::PresolveWorkingSet(solution, memspace, workingset, qrmanager);
    if (solution->state >= 0) {
      guard1 = true;
    }
  } else {
    solution->state = 82;
    guard1 = true;
  }
  if (guard1) {
    double x;
    int mLB;
    solution->iterations = 0;
    mLB = workingset->sizes[3];
    x = 0.0;
    for (idx = 0; idx < mLB; idx++) {
      idxLB = workingset->indexLB[idx] - 1;
      x = std::fmax(x, -solution->xstar[idxLB] - workingset->lb[idxLB]);
    }
    mLB = workingset->indexUB[0] - 1;
    x = std::fmax(x, solution->xstar[mLB] - workingset->ub[mLB]);
    mLB = workingset->indexUB[1] - 1;
    x = std::fmax(x, solution->xstar[mLB] - workingset->ub[mLB]);
    solution->maxConstr = x;
    if (x > 1.0E-6) {
      int PROBTYPE_ORIG;
      int b_nVar;
      int idxEndIneq;
      PROBTYPE_ORIG = workingset->probType;
      b_nVar = workingset->nVar;
      solution->xstar[4] = x + 1.0;
      if (workingset->probType == 3) {
        mLB = 1;
      } else {
        mLB = 4;
      }
      idxLB = (workingset->nWConstr[0] + workingset->nWConstr[1]) + 1;
      idxEndIneq = workingset->nActiveConstr;
      for (idx = idxLB; idx <= idxEndIneq; idx++) {
        workingset->isActiveConstr
            [(workingset->isActiveIdx[workingset->Wid[idx - 1] - 1] +
              workingset->Wlocalidx[idx - 1]) -
             2] = false;
      }
      workingset->nWConstr[2] = 0;
      workingset->nWConstr[3] = 0;
      workingset->nWConstr[4] = 0;
      workingset->nActiveConstr =
          workingset->nWConstr[0] + workingset->nWConstr[1];
      WorkingSet::setProblemType(workingset, mLB);
      objective->prev_objtype = objective->objtype;
      objective->prev_nvar = objective->nvar;
      objective->prev_hasLinear = objective->hasLinear;
      objective->objtype = 5;
      objective->nvar = 5;
      objective->gammaScalar = 1.0;
      objective->hasLinear = true;
      solution->fstar = Objective::computeFval(
          objective, memspace->workspace_double, H, f, solution->xstar);
      solution->state = 5;
      iterate(H, f, solution, memspace, workingset, qrmanager, cholmanager,
              objective, options->SolverName, 1.4901161193847657E-10, 1.0E-6,
              runTimeOptions_MaxIterations);
      if (workingset->isActiveConstr
              [(workingset->isActiveIdx[3] + workingset->sizes[3]) - 2]) {
        boolean_T exitg1;
        idx = 0;
        exitg1 = false;
        while ((!exitg1) && (idx + 1 <= workingset->nActiveConstr)) {
          if ((workingset->Wid[idx] == 4) &&
              (workingset->Wlocalidx[idx] == workingset->sizes[3])) {
            idxEndIneq = workingset->Wid[idx] - 1;
            workingset->isActiveConstr
                [(workingset->isActiveIdx[workingset->Wid[idx] - 1] +
                  workingset->Wlocalidx[idx]) -
                 2] = false;
            workingset->Wid[idx] =
                workingset->Wid[workingset->nActiveConstr - 1];
            workingset->Wlocalidx[idx] =
                workingset->Wlocalidx[workingset->nActiveConstr - 1];
            idxLB = workingset->nVar;
            for (mLB = 0; mLB < idxLB; mLB++) {
              workingset->ATwset[mLB + 5 * idx] =
                  workingset->ATwset[mLB + 5 * (workingset->nActiveConstr - 1)];
            }
            workingset->bwset[idx] =
                workingset->bwset[workingset->nActiveConstr - 1];
            workingset->nActiveConstr--;
            workingset->nWConstr[idxEndIneq]--;
            exitg1 = true;
          } else {
            idx++;
          }
        }
      }
      mLB = workingset->nActiveConstr - 1;
      while ((mLB + 1 > 0) && (mLB + 1 > b_nVar)) {
        idxEndIneq = workingset->Wid[mLB] - 1;
        workingset->isActiveConstr
            [(workingset->isActiveIdx[workingset->Wid[mLB] - 1] +
              workingset->Wlocalidx[mLB]) -
             2] = false;
        workingset->Wid[mLB] = workingset->Wid[workingset->nActiveConstr - 1];
        workingset->Wlocalidx[mLB] =
            workingset->Wlocalidx[workingset->nActiveConstr - 1];
        idxLB = workingset->nVar;
        for (idx = 0; idx < idxLB; idx++) {
          workingset->ATwset[idx + 5 * mLB] =
              workingset->ATwset[idx + 5 * (workingset->nActiveConstr - 1)];
        }
        workingset->bwset[mLB] =
            workingset->bwset[workingset->nActiveConstr - 1];
        workingset->nActiveConstr--;
        workingset->nWConstr[idxEndIneq]--;
        mLB--;
      }
      solution->maxConstr = solution->xstar[4];
      WorkingSet::setProblemType(workingset, PROBTYPE_ORIG);
      objective->objtype = objective->prev_objtype;
      objective->nvar = objective->prev_nvar;
      objective->hasLinear = objective->prev_hasLinear;
      options->ObjectiveLimit = rtMinusInf;
      options->StepTolerance = 1.0E-6;
      if (solution->state != 0) {
        mLB = workingset->sizes[3];
        x = 0.0;
        for (idx = 0; idx < mLB; idx++) {
          idxLB = workingset->indexLB[idx] - 1;
          x = std::fmax(x, -solution->xstar[idxLB] - workingset->lb[idxLB]);
        }
        mLB = workingset->indexUB[0] - 1;
        x = std::fmax(x, solution->xstar[mLB] - workingset->ub[mLB]);
        mLB = workingset->indexUB[1] - 1;
        x = std::fmax(x, solution->xstar[mLB] - workingset->ub[mLB]);
        solution->maxConstr = x;
        if (x > 1.0E-6) {
          std::memset(&solution->lambda[0], 0, 9U * sizeof(double));
          solution->fstar = Objective::computeFval(
              objective, memspace->workspace_double, H, f, solution->xstar);
          solution->state = -2;
        } else {
          if (x > 0.0) {
            if (0 <= nVar) {
              std::copy(&solution->xstar[0], &solution->xstar[nVar + 1],
                        &solution->searchDir[0]);
            }
            initialize::PresolveWorkingSet(solution, memspace, workingset,
                                           qrmanager);
            mLB = workingset->sizes[3];
            x = 0.0;
            for (idx = 0; idx < mLB; idx++) {
              idxLB = workingset->indexLB[idx] - 1;
              x = std::fmax(x, -solution->xstar[idxLB] - workingset->lb[idxLB]);
            }
            mLB = workingset->indexUB[0] - 1;
            x = std::fmax(x, solution->xstar[mLB] - workingset->ub[mLB]);
            mLB = workingset->indexUB[1] - 1;
            x = std::fmax(x, solution->xstar[mLB] - workingset->ub[mLB]);
            if (x >= solution->maxConstr) {
              solution->maxConstr = x;
              if (0 <= nVar) {
                std::copy(&solution->searchDir[0],
                          &solution->searchDir[nVar + 1], &solution->xstar[0]);
              }
            }
          }
          iterate(H, f, solution, memspace, workingset, qrmanager, cholmanager,
                  objective, options->SolverName, options->StepTolerance,
                  options->ObjectiveLimit, runTimeOptions_MaxIterations);
        }
      }
    } else {
      iterate(H, f, solution, memspace, workingset, qrmanager, cholmanager,
              objective, options->SolverName, options->StepTolerance,
              options->ObjectiveLimit, runTimeOptions_MaxIterations);
    }
  }
}

} // namespace qpactiveset
} // namespace coder
} // namespace optim
} // namespace coder

//
// File trailer for driver1.cpp
//
// [EOF]
//
