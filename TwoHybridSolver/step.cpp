//
// File: step.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

// Include Files
#include "step.h"
#include "TwoHybridSolver_internal_types.h"
#include "driver1.h"
#include "rt_nonfinite.h"
#include "setProblemType.h"
#include "sortLambdaQP.h"
#include "xnrm2.h"
#include "rtGetInf.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <string.h>

// Function Definitions
//
// Arguments    : int *STEP_TYPE
//                double Hessian[16]
//                i_struct_T *TrialState
//                b_struct_T *MeritFunction
//                d_struct_T *memspace
//                j_struct_T *WorkingSet
//                e_struct_T *b_QRManager
//                f_struct_T *b_CholManager
//                struct_T *QPObjective
//                g_struct_T *qpoptions
// Return Type  : boolean_T
//
namespace coder {
namespace optim {
namespace coder {
namespace fminconsqp {
boolean_T step(int *STEP_TYPE, double Hessian[16], i_struct_T *TrialState,
               b_struct_T *MeritFunction, d_struct_T *memspace,
               j_struct_T *WorkingSet, e_struct_T *b_QRManager,
               f_struct_T *b_CholManager, struct_T *QPObjective,
               g_struct_T *qpoptions)
{
  static double dv[4]{0.0, 0.0, 1.0, 1.0};
  g_struct_T b_qpoptions;
  double dv1[5];
  double oldDirIdx;
  double s;
  int i;
  int idxStartIneq;
  int nVar_tmp_tmp;
  boolean_T checkBoundViolation;
  boolean_T stepSuccess;
  dv[0U] = rtGetInf();
  dv[1U] = rtGetInf();
  stepSuccess = true;
  checkBoundViolation = true;
  nVar_tmp_tmp = WorkingSet->nVar - 1;
  if (*STEP_TYPE != 3) {
    if (0 <= nVar_tmp_tmp) {
      std::copy(&TrialState->xstarsqp[0],
                &TrialState->xstarsqp[nVar_tmp_tmp + 1], &TrialState->xstar[0]);
    }
  } else if (0 <= nVar_tmp_tmp) {
    std::copy(&TrialState->xstar[0], &TrialState->xstar[nVar_tmp_tmp + 1],
              &TrialState->searchDir[0]);
  }
  int exitg1;
  boolean_T guard1{false};
  do {
    int k;
    int nVarOrig;
    exitg1 = 0;
    guard1 = false;
    switch (*STEP_TYPE) {
    case 1:
      b_qpoptions = *qpoptions;
      ::coder::optim::coder::qpactiveset::driver(
          Hessian, TrialState->grad, TrialState, memspace, WorkingSet,
          b_QRManager, b_CholManager, QPObjective, &b_qpoptions,
          qpoptions->MaxIterations);
      if (TrialState->state > 0) {
        MeritFunction->phi = TrialState->sqpFval;
        MeritFunction->linearizedConstrViol = 0.0;
        MeritFunction->penaltyParam = 1.0;
        MeritFunction->phiPrimePlus = std::fmin(TrialState->fstar, 0.0);
      }
      qpactiveset::parseoutput::sortLambdaQP(
          TrialState->lambda, WorkingSet->nActiveConstr, WorkingSet->sizes,
          WorkingSet->isActiveIdx, WorkingSet->Wid, WorkingSet->Wlocalidx,
          memspace->workspace_double);
      if ((TrialState->state <= 0) && (TrialState->state != -6)) {
        *STEP_TYPE = 2;
      } else {
        if (0 <= nVar_tmp_tmp) {
          std::copy(&TrialState->xstar[0], &TrialState->xstar[nVar_tmp_tmp + 1],
                    &TrialState->delta_x[0]);
        }
        guard1 = true;
      }
      break;
    case 2: {
      double beta;
      idxStartIneq = (WorkingSet->nWConstr[0] + WorkingSet->nWConstr[1]) + 1;
      i = WorkingSet->nActiveConstr;
      for (nVarOrig = idxStartIneq; nVarOrig <= i; nVarOrig++) {
        WorkingSet->isActiveConstr
            [(WorkingSet->isActiveIdx[WorkingSet->Wid[nVarOrig - 1] - 1] +
              WorkingSet->Wlocalidx[nVarOrig - 1]) -
             2] = false;
      }
      WorkingSet->nWConstr[2] = 0;
      WorkingSet->nWConstr[3] = 0;
      WorkingSet->nWConstr[4] = 0;
      WorkingSet->nActiveConstr =
          WorkingSet->nWConstr[0] + WorkingSet->nWConstr[1];
      for (i = 0; i < 5; i++) {
        dv1[i] = TrialState->xstar[i];
      }
      idxStartIneq = WorkingSet->sizes[3];
      for (i = 0; i < idxStartIneq; i++) {
        oldDirIdx = WorkingSet->lb[WorkingSet->indexLB[i] - 1];
        if (-dv1[WorkingSet->indexLB[i] - 1] > oldDirIdx) {
          if (std::isinf(dv[WorkingSet->indexLB[i] - 1])) {
            dv1[WorkingSet->indexLB[i] - 1] = -oldDirIdx + std::abs(oldDirIdx);
          } else {
            dv1[WorkingSet->indexLB[i] - 1] =
                (WorkingSet->ub[WorkingSet->indexLB[i] - 1] - oldDirIdx) / 2.0;
          }
        }
      }
      oldDirIdx = WorkingSet->ub[WorkingSet->indexUB[0] - 1];
      if (dv1[WorkingSet->indexUB[0] - 1] > oldDirIdx) {
        dv1[WorkingSet->indexUB[0] - 1] =
            (oldDirIdx - WorkingSet->lb[WorkingSet->indexUB[0] - 1]) / 2.0;
      }
      oldDirIdx = WorkingSet->ub[WorkingSet->indexUB[1] - 1];
      if (dv1[WorkingSet->indexUB[1] - 1] > oldDirIdx) {
        dv1[WorkingSet->indexUB[1] - 1] =
            (oldDirIdx - WorkingSet->lb[WorkingSet->indexUB[1] - 1]) / 2.0;
      }
      for (i = 0; i < 5; i++) {
        TrialState->xstar[i] = dv1[i];
      }
      nVarOrig = WorkingSet->nVar;
      beta = 0.0;
      for (i = 0; i < nVarOrig; i++) {
        beta += Hessian[i + (i << 2)];
      }
      beta /= static_cast<double>(WorkingSet->nVar);
      if (TrialState->sqpIterations <= 1) {
        i = QPObjective->nvar;
        if (QPObjective->nvar < 1) {
          idxStartIneq = 0;
        } else {
          idxStartIneq = 1;
          if (QPObjective->nvar > 1) {
            oldDirIdx = std::abs(TrialState->grad[0]);
            for (k = 2; k <= i; k++) {
              s = std::abs(TrialState->grad[k - 1]);
              if (s > oldDirIdx) {
                idxStartIneq = k;
                oldDirIdx = s;
              }
            }
          }
        }
        oldDirIdx =
            100.0 *
            std::fmax(1.0, std::abs(TrialState->grad[idxStartIneq - 1]));
      } else {
        i = WorkingSet->mConstr;
        idxStartIneq = 1;
        oldDirIdx = std::abs(TrialState->lambdasqp[0]);
        for (k = 2; k <= i; k++) {
          s = std::abs(TrialState->lambdasqp[k - 1]);
          if (s > oldDirIdx) {
            idxStartIneq = k;
            oldDirIdx = s;
          }
        }
        oldDirIdx = std::abs(TrialState->lambdasqp[idxStartIneq - 1]);
      }
      QPObjective->nvar = WorkingSet->nVar;
      QPObjective->beta = beta;
      QPObjective->rho = oldDirIdx;
      QPObjective->hasLinear = true;
      QPObjective->objtype = 4;
      qpactiveset::WorkingSet::setProblemType(WorkingSet, 2);
      idxStartIneq = qpoptions->MaxIterations;
      qpoptions->MaxIterations =
          (qpoptions->MaxIterations + WorkingSet->nVar) - nVarOrig;
      for (i = 0; i < 5; i++) {
        dv1[i] = TrialState->grad[i];
      }
      b_qpoptions = *qpoptions;
      ::coder::optim::coder::qpactiveset::driver(
          Hessian, dv1, TrialState, memspace, WorkingSet, b_QRManager,
          b_CholManager, QPObjective, &b_qpoptions, qpoptions->MaxIterations);
      qpoptions->MaxIterations = idxStartIneq;
      if (TrialState->state != -6) {
        MeritFunction->phi = TrialState->sqpFval;
        MeritFunction->linearizedConstrViol = 0.0;
        MeritFunction->penaltyParam = 1.0;
        MeritFunction->phiPrimePlus = std::fmin(
            (TrialState->fstar - oldDirIdx * 0.0) - beta / 2.0 * 0.0, 0.0);
        idxStartIneq = WorkingSet->nActiveConstr;
        for (i = 1; i <= idxStartIneq; i++) {
          if (WorkingSet->Wid[i - 1] == 3) {
            TrialState->lambda[i - 1] *= static_cast<double>(
                memspace->workspace_int[WorkingSet->Wlocalidx[i - 1] - 1]);
          }
        }
      }
      QPObjective->nvar = nVarOrig;
      QPObjective->hasLinear = true;
      QPObjective->objtype = 3;
      qpactiveset::WorkingSet::setProblemType(WorkingSet, 3);
      qpactiveset::parseoutput::sortLambdaQP(
          TrialState->lambda, WorkingSet->nActiveConstr, WorkingSet->sizes,
          WorkingSet->isActiveIdx, WorkingSet->Wid, WorkingSet->Wlocalidx,
          memspace->workspace_double);
      if (0 <= nVar_tmp_tmp) {
        std::copy(&TrialState->xstar[0], &TrialState->xstar[nVar_tmp_tmp + 1],
                  &TrialState->delta_x[0]);
      }
      guard1 = true;
    } break;
    default:
      idxStartIneq = WorkingSet->nVar - 1;
      if (0 <= idxStartIneq) {
        std::copy(&TrialState->xstarsqp_old[0],
                  &TrialState->xstarsqp_old[idxStartIneq + 1],
                  &TrialState->xstarsqp[0]);
        std::copy(&TrialState->xstar[0], &TrialState->xstar[idxStartIneq + 1],
                  &TrialState->socDirection[0]);
      }
      std::copy(&TrialState->lambda[0], &TrialState->lambda[9],
                &TrialState->lambda_old[0]);
      if (0 <= idxStartIneq) {
        std::copy(&TrialState->xstarsqp[0],
                  &TrialState->xstarsqp[idxStartIneq + 1],
                  &TrialState->xstar[0]);
      }
      for (i = 0; i < 5; i++) {
        dv1[i] = TrialState->grad[i];
      }
      b_qpoptions = *qpoptions;
      ::coder::optim::coder::qpactiveset::driver(
          Hessian, dv1, TrialState, memspace, WorkingSet, b_QRManager,
          b_CholManager, QPObjective, &b_qpoptions, qpoptions->MaxIterations);
      for (i = 0; i <= idxStartIneq; i++) {
        oldDirIdx = TrialState->socDirection[i];
        TrialState->socDirection[i] =
            TrialState->xstar[i] - TrialState->socDirection[i];
        TrialState->xstar[i] = oldDirIdx;
      }
      stepSuccess = (::coder::internal::blas::xnrm2(idxStartIneq + 1,
                                                    TrialState->socDirection) <=
                     2.0 * ::coder::internal::blas::xnrm2(idxStartIneq + 1,
                                                          TrialState->xstar));
      if (!stepSuccess) {
        std::copy(&TrialState->lambda_old[0], &TrialState->lambda_old[9],
                  &TrialState->lambda[0]);
      } else {
        qpactiveset::parseoutput::sortLambdaQP(
            TrialState->lambda, WorkingSet->nActiveConstr, WorkingSet->sizes,
            WorkingSet->isActiveIdx, WorkingSet->Wid, WorkingSet->Wlocalidx,
            memspace->workspace_double);
      }
      checkBoundViolation = stepSuccess;
      if (stepSuccess && (TrialState->state != -6)) {
        for (i = 0; i <= nVar_tmp_tmp; i++) {
          TrialState->delta_x[i] =
              TrialState->xstar[i] + TrialState->socDirection[i];
        }
      }
      guard1 = true;
      break;
    }
    if (guard1) {
      if (TrialState->state != -6) {
        exitg1 = 1;
      } else {
        oldDirIdx = std::fmax(
            2.2204460492503131E-16,
            std::fmax(
                std::fmax(
                    std::fmax(std::fmax(0.0, std::abs(TrialState->grad[0])),
                              std::abs(TrialState->grad[1])),
                    std::abs(TrialState->grad[2])),
                std::abs(TrialState->grad[3])) /
                std::fmax(
                    std::fmax(
                        std::fmax(
                            std::fmax(1.0, std::abs(TrialState->xstar[0])),
                            std::abs(TrialState->xstar[1])),
                        std::abs(TrialState->xstar[2])),
                    std::abs(TrialState->xstar[3])));
        for (nVarOrig = 0; nVarOrig < 4; nVarOrig++) {
          idxStartIneq = nVarOrig << 2;
          for (k = 0; k < nVarOrig; k++) {
            Hessian[idxStartIneq + k] = 0.0;
          }
          idxStartIneq += nVarOrig;
          Hessian[idxStartIneq] = oldDirIdx;
          i = 2 - nVarOrig;
          if (0 <= i) {
            std::memset(&Hessian[idxStartIneq + 1], 0,
                        (((i + idxStartIneq) - idxStartIneq) + 1) *
                            sizeof(double));
          }
        }
      }
    }
  } while (exitg1 == 0);
  if (checkBoundViolation) {
    idxStartIneq = WorkingSet->sizes[3];
    for (i = 0; i < 5; i++) {
      dv1[i] = TrialState->delta_x[i];
    }
    for (i = 0; i < idxStartIneq; i++) {
      oldDirIdx = dv1[WorkingSet->indexLB[i] - 1];
      s = TrialState->xstarsqp[WorkingSet->indexLB[i] - 1] + oldDirIdx;
      if (s < 0.0) {
        dv1[WorkingSet->indexLB[i] - 1] = oldDirIdx - s;
        TrialState->xstar[WorkingSet->indexLB[i] - 1] -= s;
      }
    }
    oldDirIdx = dv1[WorkingSet->indexUB[0] - 1];
    s = (1.0 - TrialState->xstarsqp[WorkingSet->indexUB[0] - 1]) - oldDirIdx;
    if (s < 0.0) {
      dv1[WorkingSet->indexUB[0] - 1] = oldDirIdx + s;
      TrialState->xstar[WorkingSet->indexUB[0] - 1] += s;
    }
    oldDirIdx = dv1[WorkingSet->indexUB[1] - 1];
    s = (1.0 - TrialState->xstarsqp[WorkingSet->indexUB[1] - 1]) - oldDirIdx;
    if (s < 0.0) {
      dv1[WorkingSet->indexUB[1] - 1] = oldDirIdx + s;
      TrialState->xstar[WorkingSet->indexUB[1] - 1] += s;
    }
    for (i = 0; i < 5; i++) {
      TrialState->delta_x[i] = dv1[i];
    }
  }
  return stepSuccess;
}

} // namespace fminconsqp
} // namespace coder
} // namespace optim
} // namespace coder

//
// File trailer for step.cpp
//
// [EOF]
//
