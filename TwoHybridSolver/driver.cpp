//
// File: driver.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

// Include Files
#include "driver.h"
#include "BFGSUpdate.h"
#include "TwoHybridSolver.h"
#include "TwoHybridSolver_internal_types.h"
#include "anonymous_function.h"
#include "computeForwardDifferences.h"
#include "rt_nonfinite.h"
#include "step.h"
#include "test_exit.h"
#include "coder_array.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <string.h>

// Function Definitions
//
// Arguments    : i_struct_T *TrialState
//                b_struct_T *MeritFunction
//                const k_struct_T *FcnEvaluator
//                l_struct_T *FiniteDifferences
//                d_struct_T *memspace
//                j_struct_T *WorkingSet
//                e_struct_T *b_QRManager
//                struct_T *QPObjective
//                double Hessian[16]
//                f_struct_T *b_CholManager
// Return Type  : void
//
namespace coder {
namespace optim {
namespace coder {
namespace fminconsqp {
void driver(i_struct_T *TrialState, b_struct_T *MeritFunction,
            const k_struct_T *FcnEvaluator, l_struct_T *FiniteDifferences,
            d_struct_T *memspace, j_struct_T *WorkingSet,
            e_struct_T *b_QRManager, struct_T *QPObjective, double Hessian[16],
            f_struct_T *b_CholManager)
{
  static const signed char iv[16]{1, 0, 0, 0, 0, 1, 0, 0,
                                  0, 0, 1, 0, 0, 0, 0, 1};
  static const char qpoptions_SolverName[7]{'f', 'm', 'i', 'n', 'c', 'o', 'n'};
  c_struct_T Flags;
  g_struct_T b_expl_temp;
  g_struct_T expl_temp;
  int i;
  int mConstr;
  int mLB;
  int nVar_tmp_tmp;
  int qpoptions_MaxIterations;
  int u0;
  int u1;
  std::memset(&b_CholManager->FMat[0], 0, 81U * sizeof(double));
  b_CholManager->ldm = 9;
  b_CholManager->ndims = 0;
  b_CholManager->info = 0;
  b_CholManager->scaleFactor = 0.0;
  b_CholManager->ConvexCheck = true;
  b_CholManager->regTol_ = rtInf;
  b_CholManager->workspace_ = rtInf;
  b_CholManager->workspace2_ = rtInf;
  for (i = 0; i < 16; i++) {
    Hessian[i] = iv[i];
  }
  nVar_tmp_tmp = WorkingSet->nVar - 1;
  mLB = WorkingSet->sizes[3];
  mConstr = WorkingSet->sizes[3] + 1;
  u0 = WorkingSet->nVar;
  u1 = WorkingSet->sizes[3] + 2;
  if (u0 > u1) {
    u1 = u0;
  }
  qpoptions_MaxIterations = 10 * u1;
  TrialState->steplength = 1.0;
  test_exit(MeritFunction, WorkingSet, TrialState, &Flags.gradOK,
            &Flags.fevalOK, &Flags.done, &Flags.stepAccepted,
            &Flags.failedLineSearch, &Flags.stepType);
  TrialState->sqpFval_old = TrialState->sqpFval;
  TrialState->xstarsqp_old[0] = TrialState->xstarsqp[0];
  TrialState->grad_old[0] = TrialState->grad[0];
  TrialState->xstarsqp_old[1] = TrialState->xstarsqp[1];
  TrialState->grad_old[1] = TrialState->grad[1];
  TrialState->xstarsqp_old[2] = TrialState->xstarsqp[2];
  TrialState->grad_old[2] = TrialState->grad[2];
  TrialState->xstarsqp_old[3] = TrialState->xstarsqp[3];
  TrialState->grad_old[3] = TrialState->grad[3];
  if (!Flags.done) {
    TrialState->sqpIterations = 1;
  }
  while (!Flags.done) {
    while (!(Flags.stepAccepted || Flags.failedLineSearch)) {
      double phi_alpha;
      for (u0 = 0; u0 < mLB; u0++) {
        WorkingSet->lb[WorkingSet->indexLB[u0] - 1] =
            TrialState->xstarsqp[WorkingSet->indexLB[u0] - 1];
      }
      WorkingSet->ub[WorkingSet->indexUB[0] - 1] =
          1.0 - TrialState->xstarsqp[WorkingSet->indexUB[0] - 1];
      WorkingSet->ub[WorkingSet->indexUB[1] - 1] =
          1.0 - TrialState->xstarsqp[WorkingSet->indexUB[1] - 1];
      if (WorkingSet->nActiveConstr > 0) {
        i = WorkingSet->nActiveConstr;
        for (u0 = 1; u0 <= i; u0++) {
          switch (WorkingSet->Wid[u0 - 1]) {
          case 4:
            WorkingSet->bwset[u0 - 1] =
                WorkingSet->lb
                    [WorkingSet->indexLB[WorkingSet->Wlocalidx[u0 - 1] - 1] -
                     1];
            break;
          case 5:
            WorkingSet->bwset[u0 - 1] =
                WorkingSet->ub
                    [WorkingSet->indexUB[WorkingSet->Wlocalidx[u0 - 1] - 1] -
                     1];
            break;
          default:
            // A check that is always false is detected at compile-time.
            // Eliminating code that follows.
            break;
          }
        }
      }
      expl_temp.IterDisplayQP = false;
      expl_temp.RemainFeasible = false;
      expl_temp.ProbRelTolFactor = 1.0;
      expl_temp.ConstrRelTolFactor = 1.0;
      expl_temp.PricingTolerance = 0.0;
      expl_temp.ObjectiveLimit = rtMinusInf;
      expl_temp.ConstraintTolerance = 1.0E-6;
      expl_temp.OptimalityTolerance = 2.2204460492503131E-14;
      expl_temp.StepTolerance = 1.0E-6;
      expl_temp.MaxIterations = qpoptions_MaxIterations;
      for (i = 0; i < 7; i++) {
        expl_temp.SolverName[i] = qpoptions_SolverName[i];
      }
      b_expl_temp = expl_temp;
      Flags.stepAccepted = step(
          &Flags.stepType, Hessian, TrialState, MeritFunction, memspace,
          WorkingSet, b_QRManager, b_CholManager, QPObjective, &b_expl_temp);
      if (Flags.stepAccepted) {
        for (u0 = 0; u0 <= nVar_tmp_tmp; u0++) {
          TrialState->xstarsqp[u0] += TrialState->delta_x[u0];
        }
        phi_alpha =
            anon(FcnEvaluator->objfun.workspace.AAEST,
                 FcnEvaluator->objfun.workspace.DDEST,
                 FcnEvaluator->objfun.workspace.EACORRR,
                 FcnEvaluator->objfun.workspace.EDCORRR, TrialState->xstarsqp);
        u0 = 1;
        if (std::isinf(phi_alpha) || std::isnan(phi_alpha)) {
          if (std::isnan(phi_alpha)) {
            u0 = -3;
          } else if (phi_alpha < 0.0) {
            u0 = -1;
          } else {
            u0 = -2;
          }
        }
        TrialState->sqpFval = phi_alpha;
        if (u0 == 1) {
          u0 = 1;
        }
        Flags.fevalOK = (u0 == 1);
        TrialState->FunctionEvaluations++;
        if (Flags.fevalOK) {
          MeritFunction->phiFullStep = phi_alpha;
        } else {
          MeritFunction->phiFullStep = rtInf;
        }
      }
      if ((Flags.stepType == 1) && Flags.stepAccepted && Flags.fevalOK &&
          (MeritFunction->phi < MeritFunction->phiFullStep) &&
          (TrialState->sqpFval < TrialState->sqpFval_old)) {
        Flags.stepType = 3;
        Flags.stepAccepted = false;
      } else {
        double alpha;
        boolean_T evalWellDefined;
        boolean_T socTaken;
        if ((Flags.stepType == 3) && Flags.stepAccepted) {
          socTaken = true;
        } else {
          socTaken = false;
        }
        evalWellDefined = Flags.fevalOK;
        i = WorkingSet->nVar - 1;
        alpha = 1.0;
        u1 = 1;
        phi_alpha = MeritFunction->phiFullStep;
        if (0 <= i) {
          std::copy(&TrialState->delta_x[0], &TrialState->delta_x[i + 1],
                    &TrialState->searchDir[0]);
        }
        int exitg1;
        do {
          exitg1 = 0;
          if (TrialState->FunctionEvaluations < 400) {
            if (evalWellDefined &&
                (phi_alpha <=
                 MeritFunction->phi +
                     alpha * 0.0001 * MeritFunction->phiPrimePlus)) {
              exitg1 = 1;
            } else {
              boolean_T exitg2;
              boolean_T tooSmallX;
              alpha *= 0.7;
              for (u0 = 0; u0 <= i; u0++) {
                TrialState->delta_x[u0] = alpha * TrialState->xstar[u0];
              }
              if (socTaken) {
                phi_alpha = alpha * alpha;
                if ((i + 1 >= 1) && (!(phi_alpha == 0.0))) {
                  for (u0 = 0; u0 <= i; u0++) {
                    TrialState->delta_x[u0] +=
                        phi_alpha * TrialState->socDirection[u0];
                  }
                }
              }
              tooSmallX = true;
              u0 = 0;
              exitg2 = false;
              while ((!exitg2) && (u0 <= i)) {
                if (1.0E-6 *
                        std::fmax(1.0, std::abs(TrialState->xstarsqp[u0])) <=
                    std::abs(TrialState->delta_x[u0])) {
                  tooSmallX = false;
                  exitg2 = true;
                } else {
                  u0++;
                }
              }
              if (tooSmallX) {
                u1 = -2;
                exitg1 = 1;
              } else {
                for (u0 = 0; u0 <= i; u0++) {
                  TrialState->xstarsqp[u0] =
                      TrialState->xstarsqp_old[u0] + TrialState->delta_x[u0];
                }
                phi_alpha = anon(FcnEvaluator->objfun.workspace.AAEST,
                                 FcnEvaluator->objfun.workspace.DDEST,
                                 FcnEvaluator->objfun.workspace.EACORRR,
                                 FcnEvaluator->objfun.workspace.EDCORRR,
                                 TrialState->xstarsqp);
                u0 = 1;
                if (std::isinf(phi_alpha) || std::isnan(phi_alpha)) {
                  if (std::isnan(phi_alpha)) {
                    u0 = -3;
                  } else if (phi_alpha < 0.0) {
                    u0 = -1;
                  } else {
                    u0 = -2;
                  }
                }
                TrialState->sqpFval = phi_alpha;
                if (u0 == 1) {
                  u0 = 1;
                }
                TrialState->FunctionEvaluations++;
                evalWellDefined = (u0 == 1);
                if (!evalWellDefined) {
                  phi_alpha = rtInf;
                }
              }
            }
          } else {
            u1 = 0;
            exitg1 = 1;
          }
        } while (exitg1 == 0);
        Flags.fevalOK = evalWellDefined;
        TrialState->steplength = alpha;
        if (u1 > 0) {
          Flags.stepAccepted = true;
        } else {
          Flags.failedLineSearch = true;
        }
      }
    }
    if (Flags.stepAccepted && (!Flags.failedLineSearch)) {
      for (u0 = 0; u0 <= nVar_tmp_tmp; u0++) {
        TrialState->xstarsqp[u0] =
            TrialState->xstarsqp_old[u0] + TrialState->delta_x[u0];
      }
      for (u0 = 0; u0 <= mConstr; u0++) {
        TrialState->lambdasqp[u0] +=
            TrialState->steplength *
            (TrialState->lambda[u0] - TrialState->lambdasqp[u0]);
      }
      TrialState->sqpFval_old = TrialState->sqpFval;
      TrialState->xstarsqp_old[0] = TrialState->xstarsqp[0];
      TrialState->grad_old[0] = TrialState->grad[0];
      TrialState->xstarsqp_old[1] = TrialState->xstarsqp[1];
      TrialState->grad_old[1] = TrialState->grad[1];
      TrialState->xstarsqp_old[2] = TrialState->xstarsqp[2];
      TrialState->grad_old[2] = TrialState->grad[2];
      TrialState->xstarsqp_old[3] = TrialState->xstarsqp[3];
      TrialState->grad_old[3] = TrialState->grad[3];
      Flags.gradOK =
          utils::FiniteDifferences::internal::computeForwardDifferences(
              FiniteDifferences, TrialState->sqpFval, TrialState->xstarsqp,
              TrialState->grad);
      TrialState->FunctionEvaluations += FiniteDifferences->numEvals;
    } else {
      TrialState->sqpFval = TrialState->sqpFval_old;
      TrialState->xstarsqp[0] = TrialState->xstarsqp_old[0];
      TrialState->xstarsqp[1] = TrialState->xstarsqp_old[1];
      TrialState->xstarsqp[2] = TrialState->xstarsqp_old[2];
      TrialState->xstarsqp[3] = TrialState->xstarsqp_old[3];
    }
    b_test_exit(&Flags, memspace, MeritFunction, WorkingSet, TrialState,
                b_QRManager);
    if ((!Flags.done) && Flags.stepAccepted) {
      Flags.stepAccepted = false;
      Flags.stepType = 1;
      Flags.failedLineSearch = false;
      if (0 <= nVar_tmp_tmp) {
        std::copy(&TrialState->grad[0], &TrialState->grad[nVar_tmp_tmp + 1],
                  &TrialState->delta_gradLag[0]);
      }
      if (nVar_tmp_tmp + 1 >= 1) {
        for (u0 = 0; u0 <= nVar_tmp_tmp; u0++) {
          TrialState->delta_gradLag[u0] += -TrialState->grad_old[u0];
        }
      }
      BFGSUpdate(nVar_tmp_tmp + 1, Hessian, TrialState->delta_x,
                 TrialState->delta_gradLag, memspace->workspace_double);
      TrialState->sqpIterations++;
    }
  }
}

} // namespace fminconsqp
} // namespace coder
} // namespace optim
} // namespace coder

//
// File trailer for driver.cpp
//
// [EOF]
//
