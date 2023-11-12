//
// File: RemoveDependentIneq_.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

// Include Files
#include "RemoveDependentIneq_.h"
#include "TwoHybridSolver_internal_types.h"
#include "countsort.h"
#include "rt_nonfinite.h"
#include "xgeqp3.h"
#include <cmath>
#include <cstring>
#include <string.h>

// Function Definitions
//
// Arguments    : j_struct_T *workingset
//                e_struct_T *qrmanager
//                d_struct_T *memspace
//                double tolfactor
// Return Type  : void
//
namespace coder {
namespace optim {
namespace coder {
namespace qpactiveset {
namespace initialize {
void RemoveDependentIneq_(j_struct_T *workingset, e_struct_T *qrmanager,
                          d_struct_T *memspace, double tolfactor)
{
  int nDepIneq;
  int nFixedConstr;
  int nVar_tmp_tmp;
  nDepIneq = workingset->nActiveConstr;
  nFixedConstr = workingset->nWConstr[0] + workingset->nWConstr[1];
  nVar_tmp_tmp = workingset->nVar;
  if ((workingset->nWConstr[2] + workingset->nWConstr[3]) +
          workingset->nWConstr[4] >
      0) {
    double tol;
    int idx;
    int idxDiag;
    int idx_col;
    int iy0;
    int k;
    tol = tolfactor * static_cast<double>(workingset->nVar) *
          2.2204460492503131E-16;
    for (idx = 0; idx < nFixedConstr; idx++) {
      qrmanager->jpvt[idx] = 1;
    }
    idx_col = nFixedConstr + 1;
    if (idx_col <= nDepIneq) {
      std::memset(&qrmanager->jpvt[idx_col + -1], 0,
                  ((nDepIneq - idx_col) + 1) * sizeof(int));
    }
    for (idx_col = 0; idx_col < nDepIneq; idx_col++) {
      iy0 = 9 * idx_col;
      idxDiag = 5 * idx_col;
      for (k = 0; k < nVar_tmp_tmp; k++) {
        qrmanager->QR[iy0 + k] = workingset->ATwset[idxDiag + k];
      }
    }
    if (workingset->nVar * workingset->nActiveConstr == 0) {
      qrmanager->mrows = workingset->nVar;
      qrmanager->ncols = workingset->nActiveConstr;
      qrmanager->minRowCol = 0;
    } else {
      qrmanager->usedPivoting = true;
      qrmanager->mrows = workingset->nVar;
      qrmanager->ncols = workingset->nActiveConstr;
      idxDiag = workingset->nVar;
      iy0 = workingset->nActiveConstr;
      if (idxDiag < iy0) {
        iy0 = idxDiag;
      }
      qrmanager->minRowCol = iy0;
      internal::lapack::xgeqp3(qrmanager->QR, workingset->nVar,
                               workingset->nActiveConstr, qrmanager->jpvt,
                               qrmanager->tau);
    }
    nDepIneq = 0;
    for (idx = workingset->nActiveConstr - 1; idx + 1 > nVar_tmp_tmp; idx--) {
      nDepIneq++;
      memspace->workspace_int[nDepIneq - 1] = qrmanager->jpvt[idx];
    }
    if (idx + 1 <= workingset->nVar) {
      idxDiag = idx + 9 * idx;
      while ((idx + 1 > nFixedConstr) &&
             (std::abs(qrmanager->QR[idxDiag]) < tol)) {
        nDepIneq++;
        memspace->workspace_int[nDepIneq - 1] = qrmanager->jpvt[idx];
        idx--;
        idxDiag -= 10;
      }
    }
    utils::countsort(memspace->workspace_int, nDepIneq,
                     memspace->workspace_sort, nFixedConstr + 1,
                     workingset->nActiveConstr);
    for (idx = nDepIneq; idx >= 1; idx--) {
      iy0 = memspace->workspace_int[idx - 1] - 1;
      idxDiag = workingset->Wid[iy0] - 1;
      workingset->isActiveConstr[(workingset->isActiveIdx[idxDiag] +
                                  workingset->Wlocalidx[iy0]) -
                                 2] = false;
      workingset->Wid[iy0] = workingset->Wid[workingset->nActiveConstr - 1];
      workingset->Wlocalidx[iy0] =
          workingset->Wlocalidx[workingset->nActiveConstr - 1];
      idx_col = workingset->nVar;
      for (k = 0; k < idx_col; k++) {
        workingset->ATwset[k + 5 * iy0] =
            workingset->ATwset[k + 5 * (workingset->nActiveConstr - 1)];
      }
      workingset->bwset[iy0] = workingset->bwset[workingset->nActiveConstr - 1];
      workingset->nActiveConstr--;
      workingset->nWConstr[idxDiag]--;
    }
  }
}

} // namespace initialize
} // namespace qpactiveset
} // namespace coder
} // namespace optim
} // namespace coder

//
// File trailer for RemoveDependentIneq_.cpp
//
// [EOF]
//
