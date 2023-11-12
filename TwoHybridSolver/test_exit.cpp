//
// File: test_exit.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

// Include Files
#include "test_exit.h"
#include "TwoHybridSolver_internal_types.h"
#include "computeQ_.h"
#include "rt_nonfinite.h"
#include "sortLambdaQP.h"
#include "xgemv.h"
#include "xgeqp3.h"
#include <algorithm>
#include <cmath>
#include <string.h>

// Function Definitions
//
// Arguments    : c_struct_T *Flags
//                d_struct_T *memspace
//                b_struct_T *MeritFunction
//                const j_struct_T *WorkingSet
//                i_struct_T *TrialState
//                e_struct_T *b_QRManager
// Return Type  : void
//
namespace coder {
namespace optim {
namespace coder {
namespace fminconsqp {
void b_test_exit(c_struct_T *Flags, d_struct_T *memspace,
                 b_struct_T *MeritFunction, const j_struct_T *WorkingSet,
                 i_struct_T *TrialState, e_struct_T *b_QRManager)
{
  double feasError_tmp_tmp;
  double nlpComplErrorLSQ;
  double optimRelativeFactor;
  double s;
  double smax;
  int idx_max;
  int k;
  int mLB;
  int mLambda;
  int nVar_tmp;
  int rankR;
  boolean_T dxTooSmall;
  boolean_T exitg1;
  boolean_T isFeasible;
  nVar_tmp = WorkingSet->nVar;
  mLB = WorkingSet->sizes[3];
  mLambda = WorkingSet->sizes[3] + 1;
  idx_max = WorkingSet->sizes[3];
  if (0 <= nVar_tmp - 1) {
    std::copy(&TrialState->grad[0], &TrialState->grad[nVar_tmp],
              &TrialState->gradLag[0]);
  }
  for (rankR = 0; rankR < idx_max; rankR++) {
    TrialState->gradLag[WorkingSet->indexLB[rankR] - 1] -=
        TrialState->lambdasqp[rankR];
  }
  TrialState->gradLag[WorkingSet->indexUB[0] - 1] +=
      TrialState->lambdasqp[idx_max];
  TrialState->gradLag[WorkingSet->indexUB[1] - 1] +=
      TrialState->lambdasqp[idx_max + 1];
  if (WorkingSet->nVar < 1) {
    idx_max = 0;
  } else {
    idx_max = 1;
    if (WorkingSet->nVar > 1) {
      smax = std::abs(TrialState->grad[0]);
      for (k = 2; k <= nVar_tmp; k++) {
        s = std::abs(TrialState->grad[k - 1]);
        if (s > smax) {
          idx_max = k;
          smax = s;
        }
      }
    }
  }
  optimRelativeFactor = std::fmax(1.0, std::abs(TrialState->grad[idx_max - 1]));
  if (std::isinf(optimRelativeFactor)) {
    optimRelativeFactor = 1.0;
  }
  smax = 0.0;
  for (rankR = 0; rankR < mLB; rankR++) {
    smax = std::fmax(
        smax, 0.0 - TrialState->xstarsqp[WorkingSet->indexLB[rankR] - 1]);
  }
  s = TrialState->xstarsqp[WorkingSet->indexUB[0] - 1];
  smax = std::fmax(smax, s - 1.0);
  feasError_tmp_tmp = TrialState->xstarsqp[WorkingSet->indexUB[1] - 1];
  smax = std::fmax(smax, feasError_tmp_tmp - 1.0);
  MeritFunction->nlpPrimalFeasError = smax;
  if (TrialState->sqpIterations == 0) {
    MeritFunction->feasRelativeFactor = std::fmax(1.0, smax);
  }
  isFeasible = (smax <= 1.0E-6 * MeritFunction->feasRelativeFactor);
  dxTooSmall = true;
  nlpComplErrorLSQ = 0.0;
  rankR = 0;
  exitg1 = false;
  while ((!exitg1) && (rankR <= WorkingSet->nVar - 1)) {
    dxTooSmall = ((!std::isinf(TrialState->gradLag[rankR])) &&
                  (!std::isnan(TrialState->gradLag[rankR])));
    if (!dxTooSmall) {
      exitg1 = true;
    } else {
      nlpComplErrorLSQ =
          std::fmax(nlpComplErrorLSQ, std::abs(TrialState->gradLag[rankR]));
      rankR++;
    }
  }
  Flags->gradOK = dxTooSmall;
  MeritFunction->nlpDualFeasError = nlpComplErrorLSQ;
  if (!dxTooSmall) {
    Flags->done = true;
    if (isFeasible) {
      TrialState->sqpExitFlag = 2;
    } else {
      TrialState->sqpExitFlag = -2;
    }
  } else {
    double b_nlpComplError_tmp_tmp;
    double nlpComplErrorTmp;
    double nlpComplError_tmp_tmp;
    mLB = WorkingSet->sizes[3];
    nlpComplErrorTmp = 0.0;
    for (rankR = 0; rankR < mLB; rankR++) {
      smax = TrialState->xstarsqp[WorkingSet->indexLB[rankR] - 1];
      nlpComplErrorTmp = std::fmax(
          nlpComplErrorTmp,
          std::fmin(std::abs(smax * TrialState->lambdasqp[rankR]),
                    std::fmin(std::abs(smax), TrialState->lambdasqp[rankR])));
    }
    smax = TrialState->lambdasqp[mLB];
    nlpComplError_tmp_tmp = std::abs(1.0 - s);
    nlpComplErrorTmp = std::fmax(
        nlpComplErrorTmp, std::fmin(std::abs((1.0 - s) * smax),
                                    std::fmin(nlpComplError_tmp_tmp, smax)));
    smax = TrialState->lambdasqp[mLB + 1];
    b_nlpComplError_tmp_tmp = std::abs(1.0 - feasError_tmp_tmp);
    nlpComplErrorTmp = std::fmax(
        nlpComplErrorTmp, std::fmin(std::abs((1.0 - feasError_tmp_tmp) * smax),
                                    std::fmin(b_nlpComplError_tmp_tmp, smax)));
    MeritFunction->nlpComplError = nlpComplErrorTmp;
    s = std::fmax(nlpComplErrorLSQ, nlpComplErrorTmp);
    MeritFunction->firstOrderOpt = s;
    if (TrialState->sqpIterations > 1) {
      mLB = WorkingSet->sizes[3];
      if (0 <= nVar_tmp - 1) {
        std::copy(&TrialState->grad[0], &TrialState->grad[nVar_tmp],
                  &memspace->workspace_double[0]);
      }
      for (rankR = 0; rankR < mLB; rankR++) {
        memspace->workspace_double[WorkingSet->indexLB[rankR] - 1] -=
            TrialState->lambdasqp_old[rankR];
      }
      memspace->workspace_double[WorkingSet->indexUB[0] - 1] +=
          TrialState->lambdasqp_old[mLB];
      memspace->workspace_double[WorkingSet->indexUB[1] - 1] +=
          TrialState->lambdasqp_old[mLB + 1];
      nlpComplErrorLSQ = 0.0;
      rankR = 0;
      while ((rankR <= WorkingSet->nVar - 1) &&
             ((!std::isinf(memspace->workspace_double[rankR])) &&
              (!std::isnan(memspace->workspace_double[rankR])))) {
        nlpComplErrorLSQ = std::fmax(
            nlpComplErrorLSQ, std::abs(memspace->workspace_double[rankR]));
        rankR++;
      }
      mLB = WorkingSet->sizes[3];
      nlpComplErrorTmp = 0.0;
      for (rankR = 0; rankR < mLB; rankR++) {
        smax = TrialState->xstarsqp[WorkingSet->indexLB[rankR] - 1];
        nlpComplErrorTmp = std::fmax(
            nlpComplErrorTmp,
            std::fmin(
                std::abs(smax * TrialState->lambdasqp_old[rankR]),
                std::fmin(std::abs(smax), TrialState->lambdasqp_old[rankR])));
      }
      smax = TrialState->lambdasqp_old[mLB];
      nlpComplErrorTmp = std::fmax(
          nlpComplErrorTmp,
          std::fmin(
              std::abs(
                  (1.0 - TrialState->xstarsqp[WorkingSet->indexUB[0] - 1]) *
                  smax),
              std::fmin(nlpComplError_tmp_tmp, smax)));
      smax = TrialState->lambdasqp_old[mLB + 1];
      nlpComplErrorTmp = std::fmax(
          nlpComplErrorTmp,
          std::fmin(
              std::abs(
                  (1.0 - TrialState->xstarsqp[WorkingSet->indexUB[1] - 1]) *
                  smax),
              std::fmin(b_nlpComplError_tmp_tmp, smax)));
      smax = std::fmax(nlpComplErrorLSQ, nlpComplErrorTmp);
      if (smax < s) {
        MeritFunction->nlpDualFeasError = nlpComplErrorLSQ;
        MeritFunction->nlpComplError = nlpComplErrorTmp;
        MeritFunction->firstOrderOpt = smax;
        if (0 <= mLambda) {
          std::copy(&TrialState->lambdasqp_old[0],
                    &TrialState->lambdasqp_old[mLambda + 1],
                    &TrialState->lambdasqp[0]);
        }
      } else if (0 <= mLambda) {
        std::copy(&TrialState->lambdasqp[0],
                  &TrialState->lambdasqp[mLambda + 1],
                  &TrialState->lambdasqp_old[0]);
      }
    } else if (0 <= mLambda) {
      std::copy(&TrialState->lambdasqp[0], &TrialState->lambdasqp[mLambda + 1],
                &TrialState->lambdasqp_old[0]);
    }
    if (isFeasible &&
        (MeritFunction->nlpDualFeasError <= 1.0E-6 * optimRelativeFactor) &&
        (MeritFunction->nlpComplError <= 1.0E-6 * optimRelativeFactor)) {
      Flags->done = true;
      TrialState->sqpExitFlag = 1;
    } else {
      Flags->done = false;
      if (isFeasible && (TrialState->sqpFval < -1.0E+20)) {
        Flags->done = true;
        TrialState->sqpExitFlag = -3;
      } else {
        boolean_T guard1{false};
        guard1 = false;
        if (TrialState->sqpIterations > 0) {
          dxTooSmall = true;
          rankR = 0;
          exitg1 = false;
          while ((!exitg1) && (rankR <= nVar_tmp - 1)) {
            if (1.0E-6 *
                    std::fmax(1.0, std::abs(TrialState->xstarsqp[rankR])) <=
                std::abs(TrialState->delta_x[rankR])) {
              dxTooSmall = false;
              exitg1 = true;
            } else {
              rankR++;
            }
          }
          if (dxTooSmall) {
            if (!isFeasible) {
              if (Flags->stepType != 2) {
                Flags->stepType = 2;
                Flags->failedLineSearch = false;
                Flags->stepAccepted = false;
                guard1 = true;
              } else {
                Flags->done = true;
                TrialState->sqpExitFlag = -2;
              }
            } else {
              idx_max = WorkingSet->nActiveConstr - 1;
              if (WorkingSet->nActiveConstr > 0) {
                int minDim;
                for (k = 0; k <= idx_max; k++) {
                  TrialState->lambda[k] = 0.0;
                  mLB = 5 * k;
                  rankR = 9 * k;
                  for (minDim = 0; minDim < nVar_tmp; minDim++) {
                    b_QRManager->QR[rankR + minDim] =
                        WorkingSet->ATwset[mLB + minDim];
                  }
                }
                b_QRManager->usedPivoting = true;
                b_QRManager->mrows = WorkingSet->nVar;
                b_QRManager->ncols = WorkingSet->nActiveConstr;
                idx_max = WorkingSet->nVar;
                minDim = WorkingSet->nActiveConstr;
                if (idx_max < minDim) {
                  minDim = idx_max;
                }
                b_QRManager->minRowCol = minDim;
                ::coder::internal::lapack::xgeqp3(
                    b_QRManager->QR, WorkingSet->nVar,
                    WorkingSet->nActiveConstr, b_QRManager->jpvt,
                    b_QRManager->tau);
                QRManager::computeQ_(b_QRManager, WorkingSet->nVar);
                idx_max = WorkingSet->nVar;
                mLB = WorkingSet->nActiveConstr;
                if (idx_max > mLB) {
                  mLB = idx_max;
                }
                smax = std::abs(b_QRManager->QR[0]) *
                       std::fmin(1.4901161193847656E-8,
                                 static_cast<double>(mLB) *
                                     2.2204460492503131E-16);
                rankR = 0;
                idx_max = 0;
                while ((rankR < minDim) &&
                       (std::abs(b_QRManager->QR[idx_max]) > smax)) {
                  rankR++;
                  idx_max += 10;
                }
                ::coder::internal::blas::xgemv(
                    WorkingSet->nVar, WorkingSet->nVar, b_QRManager->Q,
                    TrialState->grad, memspace->workspace_double);
                if (rankR != 0) {
                  for (k = rankR; k >= 1; k--) {
                    idx_max = (k + (k - 1) * 9) - 1;
                    memspace->workspace_double[k - 1] /=
                        b_QRManager->QR[idx_max];
                    for (int i{0}; i <= k - 2; i++) {
                      mLB = (k - i) - 2;
                      memspace->workspace_double[mLB] -=
                          memspace->workspace_double[k - 1] *
                          b_QRManager->QR[(idx_max - i) - 1];
                    }
                  }
                }
                idx_max = WorkingSet->nActiveConstr;
                if (idx_max < minDim) {
                  minDim = idx_max;
                }
                for (rankR = 0; rankR < minDim; rankR++) {
                  TrialState->lambda[b_QRManager->jpvt[rankR] - 1] =
                      memspace->workspace_double[rankR];
                }
                qpactiveset::parseoutput::sortLambdaQP(
                    TrialState->lambda, WorkingSet->nActiveConstr,
                    WorkingSet->sizes, WorkingSet->isActiveIdx, WorkingSet->Wid,
                    WorkingSet->Wlocalidx, memspace->workspace_double);
                mLB = WorkingSet->sizes[3];
                if (0 <= nVar_tmp - 1) {
                  std::copy(&TrialState->grad[0], &TrialState->grad[nVar_tmp],
                            &memspace->workspace_double[0]);
                }
                for (rankR = 0; rankR < mLB; rankR++) {
                  memspace->workspace_double[WorkingSet->indexLB[rankR] - 1] -=
                      TrialState->lambda[rankR];
                }
                memspace->workspace_double[WorkingSet->indexUB[0] - 1] +=
                    TrialState->lambda[mLB];
                memspace->workspace_double[WorkingSet->indexUB[1] - 1] +=
                    TrialState->lambda[mLB + 1];
                s = 0.0;
                rankR = 0;
                while ((rankR <= WorkingSet->nVar - 1) &&
                       ((!std::isinf(memspace->workspace_double[rankR])) &&
                        (!std::isnan(memspace->workspace_double[rankR])))) {
                  s = std::fmax(s, std::abs(memspace->workspace_double[rankR]));
                  rankR++;
                }
                mLB = WorkingSet->sizes[3];
                nlpComplErrorLSQ = 0.0;
                for (rankR = 0; rankR < mLB; rankR++) {
                  smax = TrialState->xstarsqp[WorkingSet->indexLB[rankR] - 1];
                  nlpComplErrorLSQ = std::fmax(
                      nlpComplErrorLSQ,
                      std::fmin(std::abs(smax * TrialState->lambda[rankR]),
                                std::fmin(std::abs(smax),
                                          TrialState->lambda[rankR])));
                }
                smax = TrialState->lambda[mLB];
                nlpComplErrorLSQ = std::fmax(
                    nlpComplErrorLSQ,
                    std::fmin(
                        std::abs(
                            (1.0 -
                             TrialState->xstarsqp[WorkingSet->indexUB[0] - 1]) *
                            smax),
                        std::fmin(
                            std::abs(
                                1.0 -
                                TrialState
                                    ->xstarsqp[WorkingSet->indexUB[0] - 1]),
                            smax)));
                smax = TrialState->lambda[mLB + 1];
                nlpComplErrorLSQ = std::fmax(
                    nlpComplErrorLSQ,
                    std::fmin(
                        std::abs(
                            (1.0 -
                             TrialState->xstarsqp[WorkingSet->indexUB[1] - 1]) *
                            smax),
                        std::fmin(
                            std::abs(
                                1.0 -
                                TrialState
                                    ->xstarsqp[WorkingSet->indexUB[1] - 1]),
                            smax)));
                if ((s <= 1.0E-6 * optimRelativeFactor) &&
                    (nlpComplErrorLSQ <= 1.0E-6 * optimRelativeFactor)) {
                  MeritFunction->nlpDualFeasError = s;
                  MeritFunction->nlpComplError = nlpComplErrorLSQ;
                  MeritFunction->firstOrderOpt = std::fmax(s, nlpComplErrorLSQ);
                  if (0 <= mLambda) {
                    std::copy(&TrialState->lambda[0],
                              &TrialState->lambda[mLambda + 1],
                              &TrialState->lambdasqp[0]);
                  }
                  Flags->done = true;
                  TrialState->sqpExitFlag = 1;
                } else {
                  Flags->done = true;
                  TrialState->sqpExitFlag = 2;
                }
              } else {
                Flags->done = true;
                TrialState->sqpExitFlag = 2;
              }
            }
          } else {
            guard1 = true;
          }
        } else {
          guard1 = true;
        }
        if (guard1) {
          if (TrialState->sqpIterations >= 400) {
            Flags->done = true;
            TrialState->sqpExitFlag = 0;
          } else if (TrialState->FunctionEvaluations >= 400) {
            Flags->done = true;
            TrialState->sqpExitFlag = 0;
          }
        }
      }
    }
  }
}

//
// Arguments    : b_struct_T *MeritFunction
//                const j_struct_T *WorkingSet
//                i_struct_T *TrialState
//                boolean_T *Flags_gradOK
//                boolean_T *Flags_fevalOK
//                boolean_T *Flags_done
//                boolean_T *Flags_stepAccepted
//                boolean_T *Flags_failedLineSearch
//                int *Flags_stepType
// Return Type  : void
//
void test_exit(b_struct_T *MeritFunction, const j_struct_T *WorkingSet,
               i_struct_T *TrialState, boolean_T *Flags_gradOK,
               boolean_T *Flags_fevalOK, boolean_T *Flags_done,
               boolean_T *Flags_stepAccepted, boolean_T *Flags_failedLineSearch,
               int *Flags_stepType)
{
  double s;
  double smax;
  int idx;
  int idx_max;
  int mLB;
  int mLambda;
  int nVar_tmp;
  boolean_T exitg1;
  boolean_T isFeasible;
  *Flags_fevalOK = true;
  *Flags_done = false;
  *Flags_stepAccepted = false;
  *Flags_failedLineSearch = false;
  *Flags_stepType = 1;
  nVar_tmp = WorkingSet->nVar;
  mLB = WorkingSet->sizes[3];
  mLambda = WorkingSet->sizes[3] + 1;
  idx_max = WorkingSet->sizes[3];
  if (0 <= nVar_tmp - 1) {
    std::copy(&TrialState->grad[0], &TrialState->grad[nVar_tmp],
              &TrialState->gradLag[0]);
  }
  for (idx = 0; idx < idx_max; idx++) {
    TrialState->gradLag[WorkingSet->indexLB[idx] - 1] -=
        TrialState->lambdasqp[idx];
  }
  TrialState->gradLag[WorkingSet->indexUB[0] - 1] +=
      TrialState->lambdasqp[idx_max];
  TrialState->gradLag[WorkingSet->indexUB[1] - 1] +=
      TrialState->lambdasqp[idx_max + 1];
  if (WorkingSet->nVar < 1) {
    idx_max = 0;
  } else {
    idx_max = 1;
    if (WorkingSet->nVar > 1) {
      smax = std::abs(TrialState->grad[0]);
      for (idx = 2; idx <= nVar_tmp; idx++) {
        s = std::abs(TrialState->grad[idx - 1]);
        if (s > smax) {
          idx_max = idx;
          smax = s;
        }
      }
    }
  }
  s = std::fmax(1.0, std::abs(TrialState->grad[idx_max - 1]));
  if (std::isinf(s)) {
    s = 1.0;
  }
  smax = 0.0;
  for (idx = 0; idx < mLB; idx++) {
    smax = std::fmax(smax,
                     0.0 - TrialState->xstarsqp[WorkingSet->indexLB[idx] - 1]);
  }
  smax =
      std::fmax(smax, TrialState->xstarsqp[WorkingSet->indexUB[0] - 1] - 1.0);
  smax =
      std::fmax(smax, TrialState->xstarsqp[WorkingSet->indexUB[1] - 1] - 1.0);
  MeritFunction->nlpPrimalFeasError = smax;
  MeritFunction->feasRelativeFactor = std::fmax(1.0, smax);
  isFeasible = (smax <= 1.0E-6 * MeritFunction->feasRelativeFactor);
  *Flags_gradOK = true;
  smax = 0.0;
  idx = 0;
  exitg1 = false;
  while ((!exitg1) && (idx <= WorkingSet->nVar - 1)) {
    *Flags_gradOK = ((!std::isinf(TrialState->gradLag[idx])) &&
                     (!std::isnan(TrialState->gradLag[idx])));
    if (!*Flags_gradOK) {
      exitg1 = true;
    } else {
      smax = std::fmax(smax, std::abs(TrialState->gradLag[idx]));
      idx++;
    }
  }
  MeritFunction->nlpDualFeasError = smax;
  if (!*Flags_gradOK) {
    *Flags_done = true;
    if (isFeasible) {
      TrialState->sqpExitFlag = 2;
    } else {
      TrialState->sqpExitFlag = -2;
    }
  } else {
    MeritFunction->nlpComplError = 0.0;
    MeritFunction->firstOrderOpt = smax;
    if (0 <= mLambda) {
      std::copy(&TrialState->lambdasqp[0], &TrialState->lambdasqp[mLambda + 1],
                &TrialState->lambdasqp_old[0]);
    }
    if (isFeasible && (smax <= 1.0E-6 * s)) {
      *Flags_done = true;
      TrialState->sqpExitFlag = 1;
    } else if (isFeasible && (TrialState->sqpFval < -1.0E+20)) {
      *Flags_done = true;
      TrialState->sqpExitFlag = -3;
    } else if (TrialState->FunctionEvaluations >= 400) {
      *Flags_done = true;
      TrialState->sqpExitFlag = 0;
    }
  }
}

} // namespace fminconsqp
} // namespace coder
} // namespace optim
} // namespace coder

//
// File trailer for test_exit.cpp
//
// [EOF]
//
