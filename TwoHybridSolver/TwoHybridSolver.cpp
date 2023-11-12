//
// File: TwoHybridSolver.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

// Include Files
#include "TwoHybridSolver.h"
#include "TwoHybridSolver_internal_types.h"
#include "anonymous_function.h"
#include "computeForwardDifferences.h"
#include "driver.h"
#include "rt_nonfinite.h"
#include "setProblemType.h"
#include "sum.h"
#include "coder_array.h"
#include "rtGetInf.h"
#include <cmath>
#include <cstring>
#include <string.h>

// Function Definitions
//
// Arguments    : const coder::array<double, 2U> &AAEST
//                const coder::array<double, 2U> &DDEST
//                const coder::array<double, 2U> &EACORRR
//                const coder::array<double, 2U> &EDCORRR
//                double x[4]
//                double *y
// Return Type  : void
//
void TwoHybridSolver(const coder::array<double, 2U> &AAEST,
                     const coder::array<double, 2U> &DDEST,
                     const coder::array<double, 2U> &EACORRR,
                     const coder::array<double, 2U> &EDCORRR, double x[4],
                     double *y)
{
  static double dv[4]{0.0, 0.0, 1.0, 1.0};
  static const signed char t2_isActiveIdx[6]{1, 1, 1, 1, 5, 7};
  static const signed char t2_isActiveIdxNormal[6]{1, 1, 1, 1, 5, 7};
  static const signed char t2_isActiveIdxPhaseOne[6]{1, 1, 1, 1, 6, 8};
  static const signed char t2_isActiveIdxRegPhaseOne[6]{1, 1, 1, 1, 6, 8};
  static const signed char t2_isActiveIdxRegularized[6]{1, 1, 1, 1, 5, 7};
  static const signed char t2_indexLB[5]{1, 2, 3, 4, 0};
  static const signed char t2_indexUB[5]{3, 4, 0, 0, 0};
  static const signed char t2_sizes[5]{0, 0, 0, 4, 2};
  static const signed char t2_sizesNormal[5]{0, 0, 0, 4, 2};
  static const signed char t2_sizesPhaseOne[5]{0, 0, 0, 5, 2};
  static const signed char t2_sizesRegPhaseOne[5]{0, 0, 0, 5, 2};
  static const signed char t2_sizesRegularized[5]{0, 0, 0, 4, 2};
  b_struct_T MeritFunction;
  d_struct_T memspace;
  e_struct_T QRManager;
  f_struct_T CholManager;
  i_struct_T TrialState;
  j_struct_T WorkingSet;
  k_struct_T FcnEvaluator;
  l_struct_T FiniteDifferences;
  struct_T QPObjective;
  double Hessian[16];
  double fval;
  int b_i;
  int i;
  dv[0U] = rtGetInf();
  dv[1U] = rtGetInf();
  //  EACORRR = xlsread('test.xlsx', 'A6:A12915');
  //  EDCORRR = xlsread('test.xlsx', 'B6:B12915');
  //  AAEST = xlsread('test.xlsx', 'C6:C12915');
  //  DDEST = xlsread('test.xlsx', 'D6:D12915');
  FcnEvaluator.objfun.workspace.AAEST.set_size(1, AAEST.size(1));
  i = AAEST.size(1);
  for (b_i = 0; b_i < i; b_i++) {
    FcnEvaluator.objfun.workspace.AAEST[b_i] = AAEST[b_i];
  }
  FcnEvaluator.objfun.workspace.DDEST.set_size(1, DDEST.size(1));
  i = DDEST.size(1);
  for (b_i = 0; b_i < i; b_i++) {
    FcnEvaluator.objfun.workspace.DDEST[b_i] = DDEST[b_i];
  }
  FcnEvaluator.objfun.workspace.EACORRR.set_size(1, EACORRR.size(1));
  i = EACORRR.size(1);
  for (b_i = 0; b_i < i; b_i++) {
    FcnEvaluator.objfun.workspace.EACORRR[b_i] = EACORRR[b_i];
  }
  FcnEvaluator.objfun.workspace.EDCORRR.set_size(1, EDCORRR.size(1));
  i = EDCORRR.size(1);
  for (b_i = 0; b_i < i; b_i++) {
    FcnEvaluator.objfun.workspace.EDCORRR[b_i] = EDCORRR[b_i];
  }
  TrialState.nVarMax = 5;
  TrialState.mNonlinIneq = 0;
  TrialState.mNonlinEq = 0;
  TrialState.mIneq = 0;
  TrialState.mEq = 0;
  TrialState.iNonIneq0 = 1;
  TrialState.iNonEq0 = 1;
  TrialState.sqpFval_old = 0.0;
  TrialState.sqpIterations = 0;
  TrialState.sqpExitFlag = 0;
  std::memset(&TrialState.lambdasqp[0], 0, 9U * sizeof(double));
  TrialState.steplength = 1.0;
  for (i = 0; i < 5; i++) {
    TrialState.delta_x[i] = 0.0;
  }
  TrialState.fstar = 0.0;
  TrialState.firstorderopt = 0.0;
  std::memset(&TrialState.lambda[0], 0, 9U * sizeof(double));
  TrialState.state = 0;
  TrialState.maxConstr = 0.0;
  TrialState.iterations = 0;
  FcnEvaluator.nVar = 4;
  FcnEvaluator.mCineq = 0;
  FcnEvaluator.mCeq = 0;
  FcnEvaluator.NonFiniteSupport = true;
  FcnEvaluator.SpecifyObjectiveGradient = false;
  FcnEvaluator.SpecifyConstraintGradient = false;
  FcnEvaluator.ScaleProblem = false;
  FiniteDifferences.objfun = FcnEvaluator.objfun;
  FiniteDifferences.f_1 = 0.0;
  FiniteDifferences.f_2 = 0.0;
  FiniteDifferences.nVar = 4;
  FiniteDifferences.mIneq = 0;
  FiniteDifferences.mEq = 0;
  FiniteDifferences.numEvals = 0;
  FiniteDifferences.SpecifyObjectiveGradient = false;
  FiniteDifferences.SpecifyConstraintGradient = false;
  FiniteDifferences.isEmptyNonlcon = true;
  FiniteDifferences.FiniteDifferenceType = 0;
  FiniteDifferences.hasLB[0] = true;
  FiniteDifferences.hasUB[0] = false;
  for (i = 1; i < 4; i++) {
    FiniteDifferences.hasLB[i] = true;
    FiniteDifferences.hasUB[i] = !std::isinf(dv[i]);
  }
  FiniteDifferences.hasBounds = true;
  WorkingSet.mConstr = 6;
  WorkingSet.mConstrOrig = 6;
  WorkingSet.mConstrMax = 9;
  WorkingSet.nVar = 4;
  WorkingSet.nVarOrig = 4;
  WorkingSet.nVarMax = 5;
  WorkingSet.ldA = 5;
  for (i = 0; i < 5; i++) {
    WorkingSet.lb[i] = 0.0;
    WorkingSet.ub[i] = 0.0;
    WorkingSet.indexLB[i] = t2_indexLB[i];
    WorkingSet.indexUB[i] = t2_indexUB[i];
    WorkingSet.indexFixed[i] = 0;
  }
  WorkingSet.mEqRemoved = 0;
  std::memset(&WorkingSet.ATwset[0], 0, 45U * sizeof(double));
  WorkingSet.nActiveConstr = 0;
  std::memset(&WorkingSet.bwset[0], 0, 9U * sizeof(double));
  std::memset(&WorkingSet.maxConstrWorkspace[0], 0, 9U * sizeof(double));
  for (i = 0; i < 5; i++) {
    WorkingSet.sizes[i] = t2_sizes[i];
    WorkingSet.sizesNormal[i] = t2_sizesNormal[i];
    WorkingSet.sizesPhaseOne[i] = t2_sizesPhaseOne[i];
    WorkingSet.sizesRegularized[i] = t2_sizesRegularized[i];
    WorkingSet.sizesRegPhaseOne[i] = t2_sizesRegPhaseOne[i];
  }
  for (i = 0; i < 6; i++) {
    WorkingSet.isActiveIdx[i] = t2_isActiveIdx[i];
    WorkingSet.isActiveIdxNormal[i] = t2_isActiveIdxNormal[i];
    WorkingSet.isActiveIdxPhaseOne[i] = t2_isActiveIdxPhaseOne[i];
    WorkingSet.isActiveIdxRegularized[i] = t2_isActiveIdxRegularized[i];
    WorkingSet.isActiveIdxRegPhaseOne[i] = t2_isActiveIdxRegPhaseOne[i];
  }
  for (i = 0; i < 9; i++) {
    WorkingSet.isActiveConstr[i] = false;
    WorkingSet.Wid[i] = 0;
    WorkingSet.Wlocalidx[i] = 0;
  }
  for (i = 0; i < 5; i++) {
    WorkingSet.nWConstr[i] = 0;
  }
  WorkingSet.probType = 3;
  WorkingSet.SLACK0 = 1.0E-5;
  TrialState.xstarsqp[0] = 0.5;
  TrialState.xstarsqp[1] = 0.5;
  TrialState.xstarsqp[2] = 0.5;
  TrialState.xstarsqp[3] = 0.5;
  fval = anon(AAEST, DDEST, EACORRR, EDCORRR, TrialState.xstarsqp);
  TrialState.sqpFval = fval;
  coder::optim::coder::utils::FiniteDifferences::internal::
      computeForwardDifferences(&FiniteDifferences, fval, TrialState.xstarsqp,
                                TrialState.grad);
  TrialState.FunctionEvaluations = FiniteDifferences.numEvals + 1;
  WorkingSet.lb[WorkingSet.indexLB[0] - 1] = 0.5;
  WorkingSet.lb[WorkingSet.indexLB[1] - 1] = 0.5;
  WorkingSet.lb[WorkingSet.indexLB[2] - 1] = 0.5;
  WorkingSet.lb[WorkingSet.indexLB[3] - 1] = 0.5;
  WorkingSet.ub[WorkingSet.indexUB[0] - 1] = 0.5;
  WorkingSet.ub[WorkingSet.indexUB[1] - 1] = 0.5;
  coder::optim::coder::qpactiveset::WorkingSet::setProblemType(&WorkingSet, 3);
  for (i = 0; i < 9; i++) {
    WorkingSet.isActiveConstr[i] = false;
  }
  WorkingSet.nWConstr[0] = 0;
  WorkingSet.nWConstr[1] = 0;
  WorkingSet.nWConstr[2] = 0;
  WorkingSet.nWConstr[3] = 0;
  WorkingSet.nWConstr[4] = 0;
  WorkingSet.nActiveConstr = 0;
  MeritFunction.initFval = fval;
  MeritFunction.penaltyParam = 1.0;
  MeritFunction.threshold = 0.0001;
  MeritFunction.nPenaltyDecreases = 0;
  MeritFunction.linearizedConstrViol = 0.0;
  MeritFunction.initConstrViolationEq = 0.0;
  MeritFunction.initConstrViolationIneq = 0.0;
  MeritFunction.phi = 0.0;
  MeritFunction.phiPrimePlus = 0.0;
  MeritFunction.phiFullStep = 0.0;
  MeritFunction.feasRelativeFactor = 0.0;
  MeritFunction.nlpPrimalFeasError = 0.0;
  MeritFunction.nlpDualFeasError = 0.0;
  MeritFunction.nlpComplError = 0.0;
  MeritFunction.firstOrderOpt = 0.0;
  MeritFunction.hasObjective = true;
  QRManager.ldq = 9;
  std::memset(&QRManager.QR[0], 0, 81U * sizeof(double));
  std::memset(&QRManager.Q[0], 0, 81U * sizeof(double));
  QRManager.mrows = 0;
  QRManager.ncols = 0;
  std::memset(&QRManager.tau[0], 0, 9U * sizeof(double));
  for (i = 0; i < 9; i++) {
    QRManager.jpvt[i] = 0;
  }
  QRManager.minRowCol = 0;
  QRManager.usedPivoting = false;
  for (i = 0; i < 5; i++) {
    QPObjective.grad[i] = 0.0;
  }
  QPObjective.Hx[0] = 0.0;
  QPObjective.Hx[1] = 0.0;
  QPObjective.Hx[2] = 0.0;
  QPObjective.Hx[3] = 0.0;
  QPObjective.hasLinear = true;
  QPObjective.nvar = 4;
  QPObjective.maxVar = 5;
  QPObjective.beta = 0.0;
  QPObjective.rho = 0.0;
  QPObjective.objtype = 3;
  QPObjective.prev_objtype = 3;
  QPObjective.prev_nvar = 0;
  QPObjective.prev_hasLinear = false;
  QPObjective.gammaScalar = 0.0;
  coder::optim::coder::fminconsqp::driver(
      &TrialState, &MeritFunction, &FcnEvaluator, &FiniteDifferences, &memspace,
      &WorkingSet, &QRManager, &QPObjective, Hessian, &CholManager);
  x[0] = TrialState.xstarsqp[0];
  x[1] = TrialState.xstarsqp[1];
  x[2] = TrialState.xstarsqp[2];
  x[3] = TrialState.xstarsqp[3];
  *y = TrialState.sqpFval;
}

//
// Arguments    : const coder::array<double, 2U> &AAEST
//                const coder::array<double, 2U> &DDEST
//                const coder::array<double, 2U> &EACORRR
//                const coder::array<double, 2U> &EDCORRR
//                const double X[4]
// Return Type  : double
//
double anon(const coder::array<double, 2U> &AAEST,
            const coder::array<double, 2U> &DDEST,
            const coder::array<double, 2U> &EACORRR,
            const coder::array<double, 2U> &EDCORRR, const double X[4])
{
  coder::array<double, 2U> Afree;
  coder::array<double, 2U> Dfree;
  coder::array<double, 2U> b_y;
  coder::array<double, 2U> r;
  coder::array<double, 2U> y;
  double b_X;
  double d;
  int k;
  int nx;
  y.set_size(1, AAEST.size(1));
  nx = AAEST.size(1);
  for (k = 0; k < nx; k++) {
    y[k] = AAEST[k] * X[1];
  }
  Dfree.set_size(1, DDEST.size(1));
  nx = DDEST.size(1);
  for (k = 0; k < nx; k++) {
    Dfree[k] = DDEST[k] - X[0];
  }
  Afree.set_size(1, Dfree.size(1));
  nx = Dfree.size(1);
  for (k = 0; k < nx; k++) {
    Afree[k] = Dfree[k] - y[k];
  }
  b_y.set_size(1, Afree.size(1));
  nx = Afree.size(1);
  for (k = 0; k < nx; k++) {
    d = Afree[k];
    b_y[k] = d * d;
  }
  b_y.set_size(1, b_y.size(1));
  d = 4.0 * X[0];
  nx = b_y.size(1) - 1;
  for (k = 0; k <= nx; k++) {
    b_y[k] = b_y[k] + d * DDEST[k];
  }
  r.set_size(1, b_y.size(1));
  nx = b_y.size(1);
  for (k = 0; k < nx; k++) {
    r[k] = std::sqrt(b_y[k]);
  }
  Dfree.set_size(1, Dfree.size(1));
  nx = Dfree.size(1) - 1;
  for (k = 0; k <= nx; k++) {
    Dfree[k] = ((Dfree[k] - y[k]) + r[k]) / 2.0;
  }
  Afree.set_size(1, AAEST.size(1));
  nx = AAEST.size(1);
  for (k = 0; k < nx; k++) {
    Afree[k] = AAEST[k] - (DDEST[k] - Dfree[k]) / X[1];
  }
  Dfree.set_size(1, EACORRR.size(1));
  nx = EACORRR.size(1) - 1;
  for (k = 0; k <= nx; k++) {
    d = Dfree[k];
    d = EACORRR[k] - X[2] * d / (d + X[0]);
    Dfree[k] = d;
  }
  b_y.set_size(1, Dfree.size(1));
  nx = Dfree.size(1);
  for (k = 0; k < nx; k++) {
    d = Dfree[k];
    b_y[k] = d * d;
  }
  Afree.set_size(1, EDCORRR.size(1));
  b_X = X[0] / X[1];
  nx = EDCORRR.size(1) - 1;
  for (k = 0; k <= nx; k++) {
    d = Afree[k];
    d = EDCORRR[k] - X[3] * d / (d + b_X);
    Afree[k] = d;
  }
  y.set_size(1, Afree.size(1));
  nx = Afree.size(1);
  for (k = 0; k < nx; k++) {
    d = Afree[k];
    y[k] = d * d;
  }
  return coder::sum(b_y) + coder::sum(y);
}

//
// File trailer for TwoHybridSolver.cpp
//
// [EOF]
//
