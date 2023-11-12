//
// File: compute_deltax.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

// Include Files
#include "compute_deltax.h"
#include "TwoHybridSolver_internal_types.h"
#include "TwoHybridSolver_rtwutil.h"
#include "fullColLDL2_.h"
#include "rt_nonfinite.h"
#include "solve.h"
#include "xgemm.h"
#include "xpotrf.h"
#include <cmath>
#include <cstring>
#include <string.h>

// Function Definitions
//
// Arguments    : const double H[16]
//                i_struct_T *solution
//                d_struct_T *memspace
//                const e_struct_T *qrmanager
//                f_struct_T *cholmanager
//                const struct_T *objective
//                boolean_T alwaysPositiveDef
// Return Type  : void
//
namespace coder {
namespace optim {
namespace coder {
namespace qpactiveset {
void compute_deltax(const double H[16], i_struct_T *solution,
                    d_struct_T *memspace, const e_struct_T *qrmanager,
                    f_struct_T *cholmanager, const struct_T *objective,
                    boolean_T alwaysPositiveDef)
{
  int mNull_tmp;
  int nVar_tmp;
  nVar_tmp = qrmanager->mrows - 1;
  mNull_tmp = qrmanager->mrows - qrmanager->ncols;
  if (mNull_tmp <= 0) {
    if (0 <= nVar_tmp) {
      std::memset(&solution->searchDir[0], 0, (nVar_tmp + 1) * sizeof(double));
    }
  } else {
    int idx;
    for (idx = 0; idx <= nVar_tmp; idx++) {
      solution->searchDir[idx] = -objective->grad[idx];
    }
    if (qrmanager->ncols <= 0) {
      switch (objective->objtype) {
      case 5:
        break;
      case 3: {
        double smax;
        int idx_row;
        int ix;
        int k;
        int nVars;
        if (alwaysPositiveDef) {
          cholmanager->ndims = qrmanager->mrows;
          for (idx = 0; idx <= nVar_tmp; idx++) {
            idx_row = (nVar_tmp + 1) * idx;
            ix = 9 * idx;
            for (k = 0; k <= nVar_tmp; k++) {
              cholmanager->FMat[ix + k] = H[idx_row + k];
            }
          }
          cholmanager->info =
              internal::lapack::xpotrf(qrmanager->mrows, cholmanager->FMat);
        } else {
          cholmanager->ndims = qrmanager->mrows;
          for (idx = 0; idx <= nVar_tmp; idx++) {
            idx_row = qrmanager->mrows * idx;
            ix = 9 * idx;
            for (k = 0; k <= nVar_tmp; k++) {
              cholmanager->FMat[ix + k] = H[idx_row + k];
            }
          }
          if (qrmanager->mrows < 1) {
            nVars = -1;
          } else {
            nVars = 0;
            if (qrmanager->mrows > 1) {
              smax = std::abs(cholmanager->FMat[0]);
              for (k = 2; k <= nVar_tmp + 1; k++) {
                double s;
                s = std::abs(cholmanager->FMat[(k - 1) * 10]);
                if (s > smax) {
                  nVars = k - 1;
                  smax = s;
                }
              }
            }
          }
          cholmanager->regTol_ =
              std::fmax(std::abs(cholmanager->FMat[nVars + 9 * nVars]) *
                            2.2204460492503131E-16,
                        0.0);
          DynamicRegCholManager::fullColLDL2_(cholmanager, qrmanager->mrows);
          if (cholmanager->ConvexCheck) {
            idx = 0;
            int exitg1;
            do {
              exitg1 = 0;
              if (idx <= nVar_tmp) {
                if (cholmanager->FMat[idx + 9 * idx] <= 0.0) {
                  cholmanager->info = -idx - 1;
                  exitg1 = 1;
                } else {
                  idx++;
                }
              } else {
                cholmanager->ConvexCheck = false;
                exitg1 = 1;
              }
            } while (exitg1 == 0);
          }
        }
        if (cholmanager->info != 0) {
          solution->state = -6;
        } else if (alwaysPositiveDef) {
          CholManager::solve(cholmanager, solution->searchDir);
        } else {
          int i;
          idx_row = cholmanager->ndims - 2;
          if (cholmanager->ndims != 0) {
            for (idx = 0; idx <= idx_row + 1; idx++) {
              nVars = idx + idx * 9;
              i = idx_row - idx;
              for (k = 0; k <= i; k++) {
                ix = (idx + k) + 1;
                solution->searchDir[ix] -= solution->searchDir[idx] *
                                           cholmanager->FMat[(nVars + k) + 1];
              }
            }
          }
          i = cholmanager->ndims;
          for (idx = 0; idx < i; idx++) {
            solution->searchDir[idx] /= cholmanager->FMat[idx + 9 * idx];
          }
          idx_row = cholmanager->ndims;
          if (cholmanager->ndims != 0) {
            for (idx = idx_row; idx >= 1; idx--) {
              ix = (idx - 1) * 9;
              smax = solution->searchDir[idx - 1];
              i = idx + 1;
              for (k = idx_row; k >= i; k--) {
                smax -= cholmanager->FMat[(ix + k) - 1] *
                        solution->searchDir[k - 1];
              }
              solution->searchDir[idx - 1] = smax;
            }
          }
        }
      } break;
      default: {
        if (alwaysPositiveDef) {
          int idx_row;
          int k;
          int nVars;
          nVars = objective->nvar;
          cholmanager->ndims = objective->nvar;
          for (idx = 0; idx < nVars; idx++) {
            int ix;
            idx_row = nVars * idx;
            ix = 9 * idx;
            for (k = 0; k < nVars; k++) {
              cholmanager->FMat[ix + k] = H[idx_row + k];
            }
          }
          cholmanager->info =
              internal::lapack::xpotrf(objective->nvar, cholmanager->FMat);
          if (cholmanager->info != 0) {
            solution->state = -6;
          } else {
            double smax;
            int i;
            CholManager::solve(cholmanager, solution->searchDir);
            smax = 1.0 / objective->beta;
            idx_row = objective->nvar + 1;
            i = qrmanager->mrows;
            for (k = idx_row; k <= i; k++) {
              solution->searchDir[k - 1] *= smax;
            }
          }
        }
      } break;
      }
    } else {
      int nullStartIdx_tmp;
      nullStartIdx_tmp = 9 * qrmanager->ncols + 1;
      switch (objective->objtype) {
      case 5: {
        for (idx = 0; idx < mNull_tmp; idx++) {
          memspace->workspace_double[idx] =
              -qrmanager->Q[nVar_tmp + 9 * (qrmanager->ncols + idx)];
        }
        if (qrmanager->mrows != 0) {
          int i;
          int ix;
          if (0 <= nVar_tmp) {
            std::memset(&solution->searchDir[0], 0,
                        (nVar_tmp + 1) * sizeof(double));
          }
          ix = 0;
          i = nullStartIdx_tmp + 9 * (mNull_tmp - 1);
          for (idx = nullStartIdx_tmp; idx <= i; idx += 9) {
            int idx_row;
            idx_row = idx + nVar_tmp;
            for (int k{idx}; k <= idx_row; k++) {
              int nVars;
              nVars = k - idx;
              solution->searchDir[nVars] +=
                  qrmanager->Q[k - 1] * memspace->workspace_double[ix];
            }
            ix++;
          }
        }
      } break;
      default: {
        double smax;
        int i;
        int idx_row;
        int ix;
        int k;
        int nVars;
        switch (objective->objtype) {
        case 3:
          internal::blas::xgemm(qrmanager->mrows, mNull_tmp, qrmanager->mrows,
                                H, qrmanager->mrows, qrmanager->Q,
                                nullStartIdx_tmp, memspace->workspace_double);
          internal::blas::xgemm(mNull_tmp, mNull_tmp, qrmanager->mrows,
                                qrmanager->Q, nullStartIdx_tmp,
                                memspace->workspace_double, cholmanager->FMat);
          break;
        default:
          if (alwaysPositiveDef) {
            nVars = qrmanager->mrows;
            internal::blas::xgemm(objective->nvar, mNull_tmp, objective->nvar,
                                  H, objective->nvar, qrmanager->Q,
                                  nullStartIdx_tmp, memspace->workspace_double);
            for (ix = 0; ix < mNull_tmp; ix++) {
              i = objective->nvar + 1;
              for (idx_row = i; idx_row <= nVars; idx_row++) {
                memspace->workspace_double[(idx_row + 9 * ix) - 1] =
                    objective->beta *
                    qrmanager->Q[(idx_row + 9 * (ix + qrmanager->ncols)) - 1];
              }
            }
            internal::blas::xgemm(mNull_tmp, mNull_tmp, qrmanager->mrows,
                                  qrmanager->Q, nullStartIdx_tmp,
                                  memspace->workspace_double,
                                  cholmanager->FMat);
          }
          break;
        }
        if (alwaysPositiveDef) {
          cholmanager->ndims = mNull_tmp;
          cholmanager->info =
              internal::lapack::xpotrf(mNull_tmp, cholmanager->FMat);
        } else {
          cholmanager->ndims = mNull_tmp;
          nVars = 0;
          if (mNull_tmp > 1) {
            smax = std::abs(cholmanager->FMat[0]);
            for (k = 2; k <= mNull_tmp; k++) {
              double s;
              s = std::abs(cholmanager->FMat[(k - 1) * 10]);
              if (s > smax) {
                nVars = k - 1;
                smax = s;
              }
            }
          }
          cholmanager->regTol_ =
              std::fmax(std::abs(cholmanager->FMat[nVars + 9 * nVars]) *
                            2.2204460492503131E-16,
                        0.0);
          DynamicRegCholManager::fullColLDL2_(cholmanager, mNull_tmp);
          if (cholmanager->ConvexCheck) {
            idx = 0;
            int exitg1;
            do {
              exitg1 = 0;
              if (idx <= mNull_tmp - 1) {
                if (cholmanager->FMat[idx + 9 * idx] <= 0.0) {
                  cholmanager->info = -idx - 1;
                  exitg1 = 1;
                } else {
                  idx++;
                }
              } else {
                cholmanager->ConvexCheck = false;
                exitg1 = 1;
              }
            } while (exitg1 == 0);
          }
        }
        if (cholmanager->info != 0) {
          solution->state = -6;
        } else {
          if (qrmanager->mrows != 0) {
            if (0 <= mNull_tmp - 1) {
              std::memset(&memspace->workspace_double[0], 0,
                          mNull_tmp * sizeof(double));
            }
            i = nullStartIdx_tmp + 9 * (mNull_tmp - 1);
            for (idx = nullStartIdx_tmp; idx <= i; idx += 9) {
              smax = 0.0;
              idx_row = idx + nVar_tmp;
              for (k = idx; k <= idx_row; k++) {
                smax += qrmanager->Q[k - 1] * objective->grad[k - idx];
              }
              idx_row = div_nde_s32_floor(idx - nullStartIdx_tmp, 9);
              memspace->workspace_double[idx_row] += -smax;
            }
          }
          if (alwaysPositiveDef) {
            idx_row = cholmanager->ndims;
            if (cholmanager->ndims != 0) {
              for (idx = 0; idx < idx_row; idx++) {
                ix = idx * 9;
                smax = memspace->workspace_double[idx];
                for (k = 0; k < idx; k++) {
                  smax -=
                      cholmanager->FMat[ix + k] * memspace->workspace_double[k];
                }
                memspace->workspace_double[idx] =
                    smax / cholmanager->FMat[ix + idx];
              }
            }
            idx_row = cholmanager->ndims;
            if (cholmanager->ndims != 0) {
              for (idx = idx_row; idx >= 1; idx--) {
                nVars = (idx + (idx - 1) * 9) - 1;
                memspace->workspace_double[idx - 1] /= cholmanager->FMat[nVars];
                for (k = 0; k <= idx - 2; k++) {
                  ix = (idx - k) - 2;
                  memspace->workspace_double[ix] -=
                      memspace->workspace_double[idx - 1] *
                      cholmanager->FMat[(nVars - k) - 1];
                }
              }
            }
          } else {
            idx_row = cholmanager->ndims - 2;
            if (cholmanager->ndims != 0) {
              for (idx = 0; idx <= idx_row + 1; idx++) {
                nVars = idx + idx * 9;
                i = idx_row - idx;
                for (k = 0; k <= i; k++) {
                  ix = (idx + k) + 1;
                  memspace->workspace_double[ix] -=
                      memspace->workspace_double[idx] *
                      cholmanager->FMat[(nVars + k) + 1];
                }
              }
            }
            i = cholmanager->ndims;
            for (idx = 0; idx < i; idx++) {
              memspace->workspace_double[idx] /=
                  cholmanager->FMat[idx + 9 * idx];
            }
            idx_row = cholmanager->ndims;
            if (cholmanager->ndims != 0) {
              for (idx = idx_row; idx >= 1; idx--) {
                ix = (idx - 1) * 9;
                smax = memspace->workspace_double[idx - 1];
                i = idx + 1;
                for (k = idx_row; k >= i; k--) {
                  smax -= cholmanager->FMat[(ix + k) - 1] *
                          memspace->workspace_double[k - 1];
                }
                memspace->workspace_double[idx - 1] = smax;
              }
            }
          }
          if (qrmanager->mrows != 0) {
            if (0 <= nVar_tmp) {
              std::memset(&solution->searchDir[0], 0,
                          (nVar_tmp + 1) * sizeof(double));
            }
            ix = 0;
            i = nullStartIdx_tmp + 9 * (mNull_tmp - 1);
            for (idx = nullStartIdx_tmp; idx <= i; idx += 9) {
              idx_row = idx + nVar_tmp;
              for (k = idx; k <= idx_row; k++) {
                nVars = k - idx;
                solution->searchDir[nVars] +=
                    qrmanager->Q[k - 1] * memspace->workspace_double[ix];
              }
              ix++;
            }
          }
        }
      } break;
      }
    }
  }
}

} // namespace qpactiveset
} // namespace coder
} // namespace optim
} // namespace coder

//
// File trailer for compute_deltax.cpp
//
// [EOF]
//
