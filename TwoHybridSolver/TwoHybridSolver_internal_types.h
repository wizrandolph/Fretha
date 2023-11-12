//
// File: TwoHybridSolver_internal_types.h
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

#ifndef TWOHYBRIDSOLVER_INTERNAL_TYPES_H
#define TWOHYBRIDSOLVER_INTERNAL_TYPES_H

// Include Files
#include "TwoHybridSolver_types.h"
#include "anonymous_function.h"
#include "rtwtypes.h"

// Type Definitions
struct struct_T {
  double grad[5];
  double Hx[4];
  boolean_T hasLinear;
  int nvar;
  int maxVar;
  double beta;
  double rho;
  int objtype;
  int prev_objtype;
  int prev_nvar;
  boolean_T prev_hasLinear;
  double gammaScalar;
};

struct b_struct_T {
  double penaltyParam;
  double threshold;
  int nPenaltyDecreases;
  double linearizedConstrViol;
  double initFval;
  double initConstrViolationEq;
  double initConstrViolationIneq;
  double phi;
  double phiPrimePlus;
  double phiFullStep;
  double feasRelativeFactor;
  double nlpPrimalFeasError;
  double nlpDualFeasError;
  double nlpComplError;
  double firstOrderOpt;
  boolean_T hasObjective;
};

struct c_struct_T {
  boolean_T gradOK;
  boolean_T fevalOK;
  boolean_T done;
  boolean_T stepAccepted;
  boolean_T failedLineSearch;
  int stepType;
};

struct d_struct_T {
  double workspace_double[45];
  int workspace_int[9];
  int workspace_sort[9];
};

struct e_struct_T {
  int ldq;
  double QR[81];
  double Q[81];
  int jpvt[9];
  int mrows;
  int ncols;
  double tau[9];
  int minRowCol;
  boolean_T usedPivoting;
};

struct f_struct_T {
  double FMat[81];
  int ldm;
  int ndims;
  int info;
  double scaleFactor;
  boolean_T ConvexCheck;
  double regTol_;
  double workspace_;
  double workspace2_;
};

struct g_struct_T {
  char SolverName[7];
  int MaxIterations;
  double StepTolerance;
  double OptimalityTolerance;
  double ConstraintTolerance;
  double ObjectiveLimit;
  double PricingTolerance;
  double ConstrRelTolFactor;
  double ProbRelTolFactor;
  boolean_T RemainFeasible;
  boolean_T IterDisplayQP;
};

struct i_struct_T {
  int nVarMax;
  int mNonlinIneq;
  int mNonlinEq;
  int mIneq;
  int mEq;
  int iNonIneq0;
  int iNonEq0;
  double sqpFval;
  double sqpFval_old;
  double xstarsqp[4];
  double xstarsqp_old[4];
  double grad[5];
  double grad_old[5];
  int FunctionEvaluations;
  int sqpIterations;
  int sqpExitFlag;
  double lambdasqp[9];
  double lambdasqp_old[9];
  double steplength;
  double delta_x[5];
  double socDirection[5];
  double lambda_old[9];
  int workingset_old[9];
  double gradLag[5];
  double delta_gradLag[5];
  double xstar[5];
  double fstar;
  double firstorderopt;
  double lambda[9];
  int state;
  double maxConstr;
  int iterations;
  double searchDir[5];
};

struct j_struct_T {
  int mConstr;
  int mConstrOrig;
  int mConstrMax;
  int nVar;
  int nVarOrig;
  int nVarMax;
  int ldA;
  double lb[5];
  double ub[5];
  int indexLB[5];
  int indexUB[5];
  int indexFixed[5];
  int mEqRemoved;
  double ATwset[45];
  double bwset[9];
  int nActiveConstr;
  double maxConstrWorkspace[9];
  int sizes[5];
  int sizesNormal[5];
  int sizesPhaseOne[5];
  int sizesRegularized[5];
  int sizesRegPhaseOne[5];
  int isActiveIdx[6];
  int isActiveIdxNormal[6];
  int isActiveIdxPhaseOne[6];
  int isActiveIdxRegularized[6];
  int isActiveIdxRegPhaseOne[6];
  boolean_T isActiveConstr[9];
  int Wid[9];
  int Wlocalidx[9];
  int nWConstr[5];
  int probType;
  double SLACK0;
};

struct k_struct_T {
  coder::anonymous_function objfun;
  int nVar;
  int mCineq;
  int mCeq;
  boolean_T NonFiniteSupport;
  boolean_T SpecifyObjectiveGradient;
  boolean_T SpecifyConstraintGradient;
  boolean_T ScaleProblem;
};

struct l_struct_T {
  coder::anonymous_function objfun;
  double f_1;
  double f_2;
  int nVar;
  int mIneq;
  int mEq;
  int numEvals;
  boolean_T SpecifyObjectiveGradient;
  boolean_T SpecifyConstraintGradient;
  boolean_T isEmptyNonlcon;
  boolean_T hasLB[4];
  boolean_T hasUB[4];
  boolean_T hasBounds;
  int FiniteDifferenceType;
};

#endif
//
// File trailer for TwoHybridSolver_internal_types.h
//
// [EOF]
//
