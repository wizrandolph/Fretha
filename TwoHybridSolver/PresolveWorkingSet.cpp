//
// File: PresolveWorkingSet.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

// Include Files
#include "PresolveWorkingSet.h"
#include "RemoveDependentIneq_.h"
#include "TwoHybridSolver_internal_types.h"
#include "computeQ_.h"
#include "countsort.h"
#include "feasibleX0ForWorkingSet.h"
#include "rt_nonfinite.h"
#include "xgeqp3.h"
#include <cmath>
#include <cstring>
#include <string.h>

// Function Definitions
//
// Arguments    : i_struct_T *solution
//                d_struct_T *memspace
//                j_struct_T *workingset
//                e_struct_T *qrmanager
// Return Type  : void
//
namespace coder {
namespace optim {
namespace coder {
namespace qpactiveset {
namespace initialize {
void PresolveWorkingSet(i_struct_T *solution, d_struct_T *memspace,
                        j_struct_T *workingset, e_struct_T *qrmanager)
{
  double tol;
  int idxDiag;
  int idx_col;
  int ix;
  int ix0;
  int mTotalWorkingEq_tmp_tmp;
  int mWorkingFixed;
  int nDepInd;
  int nVar;
  solution->state = 82;
  nVar = workingset->nVar - 1;
  mWorkingFixed = workingset->nWConstr[0];
  mTotalWorkingEq_tmp_tmp = workingset->nWConstr[0] + workingset->nWConstr[1];
  nDepInd = 0;
  if (mTotalWorkingEq_tmp_tmp > 0) {
    int i;
    int k;
    int u0;
    for (ix = 0; ix < mTotalWorkingEq_tmp_tmp; ix++) {
      for (idx_col = 0; idx_col <= nVar; idx_col++) {
        qrmanager->QR[ix + 9 * idx_col] = workingset->ATwset[idx_col + 5 * ix];
      }
    }
    nDepInd = mTotalWorkingEq_tmp_tmp - workingset->nVar;
    if (0 > nDepInd) {
      nDepInd = 0;
    }
    if (0 <= nVar) {
      std::memset(&qrmanager->jpvt[0], 0, (nVar + 1) * sizeof(int));
    }
    i = mTotalWorkingEq_tmp_tmp * workingset->nVar;
    if (i == 0) {
      qrmanager->mrows = mTotalWorkingEq_tmp_tmp;
      qrmanager->ncols = workingset->nVar;
      qrmanager->minRowCol = 0;
    } else {
      qrmanager->usedPivoting = true;
      qrmanager->mrows = mTotalWorkingEq_tmp_tmp;
      qrmanager->ncols = workingset->nVar;
      idxDiag = workingset->nVar;
      if (mTotalWorkingEq_tmp_tmp < idxDiag) {
        idxDiag = mTotalWorkingEq_tmp_tmp;
      }
      qrmanager->minRowCol = idxDiag;
      internal::lapack::xgeqp3(qrmanager->QR, mTotalWorkingEq_tmp_tmp,
                               workingset->nVar, qrmanager->jpvt,
                               qrmanager->tau);
    }
    tol =
        100.0 * static_cast<double>(workingset->nVar) * 2.2204460492503131E-16;
    u0 = workingset->nVar;
    if (u0 >= mTotalWorkingEq_tmp_tmp) {
      u0 = mTotalWorkingEq_tmp_tmp;
    }
    idxDiag = u0 + 9 * (u0 - 1);
    while ((idxDiag > 0) && (std::abs(qrmanager->QR[idxDiag - 1]) < tol)) {
      idxDiag -= 10;
      nDepInd++;
    }
    if (nDepInd > 0) {
      boolean_T exitg1;
      QRManager::computeQ_(qrmanager, qrmanager->mrows);
      idx_col = 0;
      exitg1 = false;
      while ((!exitg1) && (idx_col <= nDepInd - 1)) {
        double qtb;
        ix = 9 * ((mTotalWorkingEq_tmp_tmp - idx_col) - 1);
        qtb = 0.0;
        for (k = 0; k < mTotalWorkingEq_tmp_tmp; k++) {
          qtb += qrmanager->Q[ix + k] * workingset->bwset[k];
        }
        if (std::abs(qtb) >= tol) {
          nDepInd = -1;
          exitg1 = true;
        } else {
          idx_col++;
        }
      }
    }
    if (nDepInd > 0) {
      for (idx_col = 0; idx_col < mTotalWorkingEq_tmp_tmp; idx_col++) {
        idxDiag = 9 * idx_col;
        ix0 = 5 * idx_col;
        for (k = 0; k <= nVar; k++) {
          qrmanager->QR[idxDiag + k] = workingset->ATwset[ix0 + k];
        }
      }
      for (idx_col = 0; idx_col < mWorkingFixed; idx_col++) {
        qrmanager->jpvt[idx_col] = 1;
      }
      ix0 = workingset->nWConstr[0] + 1;
      if (ix0 <= mTotalWorkingEq_tmp_tmp) {
        std::memset(&qrmanager->jpvt[ix0 + -1], 0,
                    ((mTotalWorkingEq_tmp_tmp - ix0) + 1) * sizeof(int));
      }
      if (i == 0) {
        qrmanager->mrows = workingset->nVar;
        qrmanager->ncols = mTotalWorkingEq_tmp_tmp;
        qrmanager->minRowCol = 0;
      } else {
        qrmanager->usedPivoting = true;
        qrmanager->mrows = workingset->nVar;
        qrmanager->ncols = mTotalWorkingEq_tmp_tmp;
        qrmanager->minRowCol = u0;
        internal::lapack::xgeqp3(qrmanager->QR, workingset->nVar,
                                 mTotalWorkingEq_tmp_tmp, qrmanager->jpvt,
                                 qrmanager->tau);
      }
      for (idx_col = 0; idx_col < nDepInd; idx_col++) {
        memspace->workspace_int[idx_col] =
            qrmanager->jpvt[(mTotalWorkingEq_tmp_tmp - nDepInd) + idx_col];
      }
      utils::countsort(memspace->workspace_int, nDepInd,
                       memspace->workspace_sort, 1, mTotalWorkingEq_tmp_tmp);
      for (idx_col = nDepInd; idx_col >= 1; idx_col--) {
        i = workingset->nWConstr[0] + workingset->nWConstr[1];
        if (i != 0) {
          ix0 = memspace->workspace_int[idx_col - 1];
          if (ix0 <= i) {
            if ((workingset->nActiveConstr == i) || (ix0 == i)) {
              workingset->mEqRemoved++;
              i = memspace->workspace_int[idx_col - 1] - 1;
              idxDiag = workingset->Wid[i] - 1;
              workingset->isActiveConstr[(workingset->isActiveIdx[idxDiag] +
                                          workingset->Wlocalidx[i]) -
                                         2] = false;
              workingset->Wid[i] =
                  workingset->Wid[workingset->nActiveConstr - 1];
              workingset->Wlocalidx[i] =
                  workingset->Wlocalidx[workingset->nActiveConstr - 1];
              ix0 = workingset->nVar;
              for (u0 = 0; u0 < ix0; u0++) {
                workingset
                    ->ATwset[u0 +
                             5 * (memspace->workspace_int[idx_col - 1] - 1)] =
                    workingset
                        ->ATwset[u0 + 5 * (workingset->nActiveConstr - 1)];
              }
              workingset->bwset[i] =
                  workingset->bwset[workingset->nActiveConstr - 1];
              workingset->nActiveConstr--;
              workingset->nWConstr[idxDiag]--;
            } else {
              workingset->mEqRemoved++;
              ix = workingset->Wid[ix0 - 1] - 1;
              workingset->isActiveConstr[(workingset->isActiveIdx[ix] +
                                          workingset->Wlocalidx[ix0 - 1]) -
                                         2] = false;
              workingset->Wid[ix0 - 1] = workingset->Wid[i - 1];
              workingset->Wlocalidx[ix0 - 1] = workingset->Wlocalidx[i - 1];
              idxDiag = workingset->nVar;
              for (u0 = 0; u0 < idxDiag; u0++) {
                workingset->ATwset[u0 + 5 * (ix0 - 1)] =
                    workingset->ATwset[u0 + 5 * (i - 1)];
              }
              workingset->bwset[ix0 - 1] = workingset->bwset[i - 1];
              workingset->Wid[i - 1] =
                  workingset->Wid[workingset->nActiveConstr - 1];
              workingset->Wlocalidx[i - 1] =
                  workingset->Wlocalidx[workingset->nActiveConstr - 1];
              ix0 = workingset->nVar;
              for (u0 = 0; u0 < ix0; u0++) {
                workingset->ATwset[u0 + 5 * (i - 1)] =
                    workingset
                        ->ATwset[u0 + 5 * (workingset->nActiveConstr - 1)];
              }
              workingset->bwset[i - 1] =
                  workingset->bwset[workingset->nActiveConstr - 1];
              workingset->nActiveConstr--;
              workingset->nWConstr[ix]--;
            }
          }
        }
      }
    }
  }
  if ((nDepInd != -1) && (workingset->nActiveConstr <= 9)) {
    boolean_T guard1{false};
    boolean_T okWorkingSet;
    RemoveDependentIneq_(workingset, qrmanager, memspace, 100.0);
    okWorkingSet = feasibleX0ForWorkingSet(
        memspace->workspace_double, solution->xstar, workingset, qrmanager);
    guard1 = false;
    if (!okWorkingSet) {
      RemoveDependentIneq_(workingset, qrmanager, memspace, 1000.0);
      okWorkingSet = feasibleX0ForWorkingSet(
          memspace->workspace_double, solution->xstar, workingset, qrmanager);
      if (!okWorkingSet) {
        solution->state = -7;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }
    if (guard1 && (workingset->nWConstr[0] + workingset->nWConstr[1] ==
                   workingset->nVar)) {
      idxDiag = workingset->sizes[3];
      tol = 0.0;
      for (idx_col = 0; idx_col < idxDiag; idx_col++) {
        ix = workingset->indexLB[idx_col] - 1;
        tol = std::fmax(tol, -solution->xstar[ix] - workingset->lb[ix]);
      }
      idxDiag = workingset->indexUB[0] - 1;
      tol = std::fmax(tol, solution->xstar[idxDiag] - workingset->ub[idxDiag]);
      idxDiag = workingset->indexUB[1] - 1;
      tol = std::fmax(tol, solution->xstar[idxDiag] - workingset->ub[idxDiag]);
      if (tol > 1.0E-6) {
        solution->state = -2;
      }
    }
  } else {
    solution->state = -3;
    idxDiag = (workingset->nWConstr[0] + workingset->nWConstr[1]) + 1;
    ix = workingset->nActiveConstr;
    for (ix0 = idxDiag; ix0 <= ix; ix0++) {
      workingset->isActiveConstr
          [(workingset->isActiveIdx[workingset->Wid[ix0 - 1] - 1] +
            workingset->Wlocalidx[ix0 - 1]) -
           2] = false;
    }
    workingset->nWConstr[2] = 0;
    workingset->nWConstr[3] = 0;
    workingset->nWConstr[4] = 0;
    workingset->nActiveConstr =
        workingset->nWConstr[0] + workingset->nWConstr[1];
  }
}

} // namespace initialize
} // namespace qpactiveset
} // namespace coder
} // namespace optim
} // namespace coder

//
// File trailer for PresolveWorkingSet.cpp
//
// [EOF]
//
