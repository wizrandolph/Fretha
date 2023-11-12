//
// File: feasibleX0ForWorkingSet.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

// Include Files
#include "feasibleX0ForWorkingSet.h"
#include "TwoHybridSolver_internal_types.h"
#include "TwoHybridSolver_rtwutil.h"
#include "computeQ_.h"
#include "factorQR.h"
#include "rt_nonfinite.h"
#include "xzgeqp3.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <string.h>

// Function Definitions
//
// Arguments    : double workspace[45]
//                double xCurrent[5]
//                const j_struct_T *workingset
//                e_struct_T *qrmanager
// Return Type  : boolean_T
//
namespace coder {
namespace optim {
namespace coder {
namespace qpactiveset {
namespace initialize {
boolean_T feasibleX0ForWorkingSet(double workspace[45], double xCurrent[5],
                                  const j_struct_T *workingset,
                                  e_struct_T *qrmanager)
{
  double B[45];
  int mWConstr;
  int nVar;
  boolean_T nonDegenerateWset;
  mWConstr = workingset->nActiveConstr;
  nVar = workingset->nVar;
  nonDegenerateWset = true;
  if (mWConstr != 0) {
    double c;
    int i;
    int i1;
    int iAcol;
    int k;
    int kAcol;
    for (kAcol = 0; kAcol < mWConstr; kAcol++) {
      workspace[kAcol] = workingset->bwset[kAcol];
      workspace[kAcol + 9] = workingset->bwset[kAcol];
    }
    if (mWConstr != 0) {
      i = 5 * (mWConstr - 1) + 1;
      for (iAcol = 1; iAcol <= i; iAcol += 5) {
        c = 0.0;
        i1 = (iAcol + nVar) - 1;
        for (kAcol = iAcol; kAcol <= i1; kAcol++) {
          c += workingset->ATwset[kAcol - 1] * xCurrent[kAcol - iAcol];
        }
        i1 = div_nde_s32_floor(iAcol - 1, 5);
        workspace[i1] += -c;
      }
    }
    if (mWConstr >= nVar) {
      int ar;
      int jBcol;
      qrmanager->usedPivoting = false;
      qrmanager->mrows = mWConstr;
      qrmanager->ncols = nVar;
      for (kAcol = 0; kAcol < nVar; kAcol++) {
        iAcol = 9 * kAcol;
        for (jBcol = 0; jBcol < mWConstr; jBcol++) {
          qrmanager->QR[jBcol + iAcol] = workingset->ATwset[kAcol + 5 * jBcol];
        }
        qrmanager->jpvt[kAcol] = kAcol + 1;
      }
      if (mWConstr < nVar) {
        i = mWConstr;
      } else {
        i = nVar;
      }
      qrmanager->minRowCol = i;
      std::memset(&qrmanager->tau[0], 0, 9U * sizeof(double));
      if (i >= 1) {
        std::memset(&qrmanager->tau[0], 0, 9U * sizeof(double));
        internal::reflapack::qrf(qrmanager->QR, mWConstr, nVar, i,
                                 qrmanager->tau);
      }
      QRManager::computeQ_(qrmanager, mWConstr);
      std::copy(&workspace[0], &workspace[45], &B[0]);
      for (k = 0; k <= 9; k += 9) {
        i = k + 1;
        i1 = k + nVar;
        if (i <= i1) {
          std::memset(&workspace[i + -1], 0, ((i1 - i) + 1) * sizeof(double));
        }
      }
      jBcol = -1;
      for (k = 0; k <= 9; k += 9) {
        ar = -1;
        i = k + 1;
        i1 = k + nVar;
        for (int ic{i}; ic <= i1; ic++) {
          c = 0.0;
          for (iAcol = 0; iAcol < mWConstr; iAcol++) {
            c += qrmanager->Q[(iAcol + ar) + 1] * B[(iAcol + jBcol) + 1];
          }
          workspace[ic - 1] += c;
          ar += 9;
        }
        jBcol += 9;
      }
      for (ar = 0; ar < 2; ar++) {
        jBcol = 9 * ar - 1;
        for (k = nVar; k >= 1; k--) {
          kAcol = 9 * (k - 1) - 1;
          i = k + jBcol;
          c = workspace[i];
          if (c != 0.0) {
            workspace[i] = c / qrmanager->QR[k + kAcol];
            for (int b_i{0}; b_i <= k - 2; b_i++) {
              i1 = (b_i + jBcol) + 1;
              workspace[i1] -= workspace[i] * qrmanager->QR[(b_i + kAcol) + 1];
            }
          }
        }
      }
    } else {
      int ar;
      int b_i;
      int jBcol;
      QRManager::factorQR(qrmanager, workingset->ATwset, nVar, mWConstr);
      QRManager::computeQ_(qrmanager, qrmanager->minRowCol);
      for (ar = 0; ar < 2; ar++) {
        jBcol = 9 * ar;
        for (b_i = 0; b_i < mWConstr; b_i++) {
          iAcol = 9 * b_i;
          kAcol = b_i + jBcol;
          c = workspace[kAcol];
          for (k = 0; k < b_i; k++) {
            c -= qrmanager->QR[k + iAcol] * workspace[k + jBcol];
          }
          workspace[kAcol] = c / qrmanager->QR[b_i + iAcol];
        }
      }
      std::copy(&workspace[0], &workspace[45], &B[0]);
      for (k = 0; k <= 9; k += 9) {
        i = k + 1;
        i1 = k + nVar;
        if (i <= i1) {
          std::memset(&workspace[i + -1], 0, ((i1 - i) + 1) * sizeof(double));
        }
      }
      jBcol = 0;
      for (k = 0; k <= 9; k += 9) {
        ar = -1;
        i = jBcol + 1;
        i1 = jBcol + mWConstr;
        for (b_i = i; b_i <= i1; b_i++) {
          iAcol = k + 1;
          kAcol = k + nVar;
          for (int ic{iAcol}; ic <= kAcol; ic++) {
            workspace[ic - 1] += B[b_i - 1] * qrmanager->Q[(ar + ic) - k];
          }
          ar += 9;
        }
        jBcol += 9;
      }
    }
    kAcol = 0;
    int exitg1;
    do {
      exitg1 = 0;
      if (kAcol <= nVar - 1) {
        if (std::isinf(workspace[kAcol]) || std::isnan(workspace[kAcol])) {
          nonDegenerateWset = false;
          exitg1 = 1;
        } else {
          c = workspace[kAcol + 9];
          if (std::isinf(c) || std::isnan(c)) {
            nonDegenerateWset = false;
            exitg1 = 1;
          } else {
            kAcol++;
          }
        }
      } else {
        double constrViolation_basicX;
        iAcol = nVar - 1;
        for (k = 0; k <= iAcol; k++) {
          workspace[k] += xCurrent[k];
        }
        iAcol = workingset->sizes[3];
        c = 0.0;
        for (kAcol = 0; kAcol < iAcol; kAcol++) {
          c = std::fmax(c, -workspace[workingset->indexLB[kAcol] - 1] -
                               workingset->lb[workingset->indexLB[kAcol] - 1]);
        }
        iAcol = workingset->indexUB[0] - 1;
        c = std::fmax(c, workspace[iAcol] - workingset->ub[iAcol]);
        iAcol = workingset->indexUB[1] - 1;
        c = std::fmax(c, workspace[iAcol] - workingset->ub[iAcol]);
        iAcol = workingset->sizes[3];
        constrViolation_basicX = 0.0;
        for (kAcol = 0; kAcol < iAcol; kAcol++) {
          constrViolation_basicX =
              std::fmax(constrViolation_basicX,
                        -workspace[workingset->indexLB[kAcol] + 8] -
                            workingset->lb[workingset->indexLB[kAcol] - 1]);
        }
        iAcol = workingset->indexUB[0];
        constrViolation_basicX =
            std::fmax(constrViolation_basicX,
                      workspace[iAcol + 8] - workingset->ub[iAcol - 1]);
        iAcol = workingset->indexUB[1];
        constrViolation_basicX =
            std::fmax(constrViolation_basicX,
                      workspace[iAcol + 8] - workingset->ub[iAcol - 1]);
        if ((c <= 2.2204460492503131E-16) || (c < constrViolation_basicX)) {
          if (0 <= nVar - 1) {
            std::copy(&workspace[0], &workspace[nVar], &xCurrent[0]);
          }
        } else if (0 <= nVar - 1) {
          std::copy(&workspace[9], &workspace[9 + nVar], &xCurrent[0]);
        }
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }
  return nonDegenerateWset;
}

} // namespace initialize
} // namespace qpactiveset
} // namespace coder
} // namespace optim
} // namespace coder

//
// File trailer for feasibleX0ForWorkingSet.cpp
//
// [EOF]
//
