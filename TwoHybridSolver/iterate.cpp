//
// File: iterate.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

// Include Files
#include "iterate.h"
#include "TwoHybridSolver_internal_types.h"
#include "TwoHybridSolver_rtwutil.h"
#include "addBoundToActiveSetMatrix_.h"
#include "checkStoppingAndUpdateFval.h"
#include "computeFval_ReuseHx.h"
#include "computeGrad_StoreHx.h"
#include "computeQ_.h"
#include "compute_deltax.h"
#include "deleteColMoveEnd.h"
#include "factorQR.h"
#include "feasibleratiotest.h"
#include "rt_nonfinite.h"
#include "xgemv.h"
#include "xnrm2.h"
#include "xrotg.h"
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
//                const char options_SolverName[7]
//                double options_StepTolerance
//                double options_ObjectiveLimit
//                int runTimeOptions_MaxIterations
// Return Type  : void
//
namespace coder {
namespace optim {
namespace coder {
namespace qpactiveset {
void iterate(const double H[16], const double f[5], i_struct_T *solution,
             d_struct_T *memspace, j_struct_T *workingset,
             e_struct_T *qrmanager, f_struct_T *cholmanager,
             struct_T *objective, const char options_SolverName[7],
             double options_StepTolerance, double options_ObjectiveLimit,
             int runTimeOptions_MaxIterations)
{
  static const char b[7]{'f', 'm', 'i', 'n', 'c', 'o', 'n'};
  double c;
  double minLambda;
  double s;
  double temp_tmp;
  int TYPE;
  int activeSetChangeID;
  int globalActiveConstrIdx;
  int iyend;
  int n;
  int nVar_tmp_tmp;
  boolean_T subProblemChanged;
  boolean_T updateFval;
  subProblemChanged = true;
  updateFval = true;
  activeSetChangeID = 0;
  TYPE = objective->objtype;
  nVar_tmp_tmp = workingset->nVar;
  globalActiveConstrIdx = 0;
  Objective::computeGrad_StoreHx(objective, H, f, solution->xstar);
  solution->fstar = Objective::computeFval_ReuseHx(
      objective, memspace->workspace_double, f, solution->xstar);
  if (solution->iterations < runTimeOptions_MaxIterations) {
    solution->state = -5;
  } else {
    solution->state = 0;
  }
  std::memset(&solution->lambda[0], 0, 9U * sizeof(double));
  int exitg1;
  do {
    exitg1 = 0;
    if (solution->state == -5) {
      int ia;
      int idx;
      int idxRotGCol;
      int nActiveConstr;
      boolean_T guard1{false};
      guard1 = false;
      if (subProblemChanged) {
        switch (activeSetChangeID) {
        case 1:
          nActiveConstr = 5 * (workingset->nActiveConstr - 1);
          iyend = qrmanager->mrows;
          idxRotGCol = qrmanager->ncols + 1;
          if (iyend < idxRotGCol) {
            idxRotGCol = iyend;
          }
          qrmanager->minRowCol = idxRotGCol;
          idxRotGCol = 9 * qrmanager->ncols;
          if (qrmanager->mrows != 0) {
            iyend = idxRotGCol + qrmanager->mrows;
            if (idxRotGCol + 1 <= iyend) {
              std::memset(&qrmanager->QR[idxRotGCol], 0,
                          (iyend - idxRotGCol) * sizeof(double));
            }
            n = 9 * (qrmanager->mrows - 1) + 1;
            for (idx = 1; idx <= n; idx += 9) {
              c = 0.0;
              iyend = (idx + qrmanager->mrows) - 1;
              for (ia = idx; ia <= iyend; ia++) {
                c += qrmanager->Q[ia - 1] *
                     workingset->ATwset[(nActiveConstr + ia) - idx];
              }
              iyend = idxRotGCol + div_nde_s32_floor(idx - 1, 9);
              qrmanager->QR[iyend] += c;
            }
          }
          qrmanager->ncols++;
          qrmanager->jpvt[qrmanager->ncols - 1] = qrmanager->ncols;
          for (idx = qrmanager->mrows - 2; idx + 2 > qrmanager->ncols; idx--) {
            idxRotGCol = 9 * (qrmanager->ncols - 1);
            n = (idx + idxRotGCol) + 1;
            temp_tmp = qrmanager->QR[n];
            internal::blas::xrotg(&qrmanager->QR[idx + idxRotGCol], &temp_tmp,
                                  &c, &s);
            qrmanager->QR[n] = temp_tmp;
            iyend = 9 * idx;
            n = qrmanager->mrows;
            if (qrmanager->mrows >= 1) {
              for (nActiveConstr = 0; nActiveConstr < n; nActiveConstr++) {
                idxRotGCol = iyend + nActiveConstr;
                minLambda = qrmanager->Q[idxRotGCol + 9];
                temp_tmp = qrmanager->Q[idxRotGCol];
                qrmanager->Q[idxRotGCol + 9] = c * minLambda - s * temp_tmp;
                qrmanager->Q[idxRotGCol] = c * temp_tmp + s * minLambda;
              }
            }
          }
          break;
        case -1:
          QRManager::deleteColMoveEnd(qrmanager, globalActiveConstrIdx);
          break;
        default:
          QRManager::factorQR(qrmanager, workingset->ATwset, nVar_tmp_tmp,
                              workingset->nActiveConstr);
          QRManager::computeQ_(qrmanager, qrmanager->mrows);
          break;
        }
        iyend = memcmp(&options_SolverName[0], &b[0], 7);
        compute_deltax(H, solution, memspace, qrmanager, cholmanager, objective,
                       iyend == 0);
        if (solution->state != -5) {
          exitg1 = 1;
        } else if ((internal::blas::xnrm2(nVar_tmp_tmp, solution->searchDir) <
                    options_StepTolerance) ||
                   (workingset->nActiveConstr >= nVar_tmp_tmp)) {
          guard1 = true;
        } else {
          feasibleratiotest(solution->xstar, solution->searchDir,
                            workingset->nVar, workingset->lb, workingset->ub,
                            workingset->indexLB, workingset->indexUB,
                            workingset->sizes, workingset->isActiveIdx,
                            workingset->isActiveConstr, workingset->nWConstr,
                            TYPE == 5, &minLambda, &updateFval, &n, &iyend);
          if (updateFval) {
            switch (n) {
            case 3:
              workingset->nWConstr[2]++;
              workingset
                  ->isActiveConstr[(workingset->isActiveIdx[2] + iyend) - 2] =
                  true;
              workingset->nActiveConstr++;
              workingset->Wid[workingset->nActiveConstr - 1] = 3;
              workingset->Wlocalidx[workingset->nActiveConstr - 1] = iyend;
              // A check that is always false is detected at compile-time.
              // Eliminating code that follows.
              break;
            case 4:
              WorkingSet::addBoundToActiveSetMatrix_(workingset, 4, iyend);
              break;
            default:
              WorkingSet::addBoundToActiveSetMatrix_(workingset, 5, iyend);
              break;
            }
            activeSetChangeID = 1;
          } else {
            if (objective->objtype == 5) {
              if (internal::blas::xnrm2(objective->nvar, solution->searchDir) >
                  100.0 * static_cast<double>(objective->nvar) *
                      1.4901161193847656E-8) {
                solution->state = 3;
              } else {
                solution->state = 4;
              }
            }
            subProblemChanged = false;
            if (workingset->nActiveConstr == 0) {
              solution->state = 1;
            }
          }
          if ((nVar_tmp_tmp >= 1) && (!(minLambda == 0.0))) {
            iyend = nVar_tmp_tmp - 1;
            for (nActiveConstr = 0; nActiveConstr <= iyend; nActiveConstr++) {
              solution->xstar[nActiveConstr] +=
                  minLambda * solution->searchDir[nActiveConstr];
            }
          }
          Objective::computeGrad_StoreHx(objective, H, f, solution->xstar);
          updateFval = true;
          stopping::checkStoppingAndUpdateFval(
              &activeSetChangeID, f, solution, memspace, objective, workingset,
              qrmanager, options_ObjectiveLimit, runTimeOptions_MaxIterations,
              updateFval);
        }
      } else {
        if (0 <= nVar_tmp_tmp - 1) {
          std::memset(&solution->searchDir[0], 0,
                      nVar_tmp_tmp * sizeof(double));
        }
        guard1 = true;
      }
      if (guard1) {
        nActiveConstr = qrmanager->ncols;
        if (qrmanager->ncols > 0) {
          minLambda = 100.0 * static_cast<double>(qrmanager->mrows) *
                      2.2204460492503131E-16;
          if ((qrmanager->mrows > 0) && (qrmanager->ncols > 0)) {
            updateFval = true;
          } else {
            updateFval = false;
          }
          if (updateFval) {
            boolean_T b_guard1{false};
            idx = qrmanager->ncols;
            b_guard1 = false;
            if (qrmanager->mrows < qrmanager->ncols) {
              iyend = qrmanager->mrows + 9 * (qrmanager->ncols - 1);
              while ((idx > qrmanager->mrows) &&
                     (std::abs(qrmanager->QR[iyend - 1]) >= minLambda)) {
                idx--;
                iyend -= 9;
              }
              updateFval = (idx == qrmanager->mrows);
              if (updateFval) {
                b_guard1 = true;
              }
            } else {
              b_guard1 = true;
            }
            if (b_guard1) {
              iyend = idx + 9 * (idx - 1);
              while ((idx >= 1) &&
                     (std::abs(qrmanager->QR[iyend - 1]) >= minLambda)) {
                idx--;
                iyend -= 10;
              }
              updateFval = (idx == 0);
            }
          }
          if (!updateFval) {
            solution->state = -7;
          } else {
            n = qrmanager->ncols;
            internal::blas::xgemv(qrmanager->mrows, qrmanager->ncols,
                                  qrmanager->Q, objective->grad,
                                  memspace->workspace_double);
            if (qrmanager->ncols != 0) {
              for (idx = n; idx >= 1; idx--) {
                iyend = (idx + (idx - 1) * 9) - 1;
                memspace->workspace_double[idx - 1] /= qrmanager->QR[iyend];
                for (ia = 0; ia <= idx - 2; ia++) {
                  idxRotGCol = (idx - ia) - 2;
                  memspace->workspace_double[idxRotGCol] -=
                      memspace->workspace_double[idx - 1] *
                      qrmanager->QR[(iyend - ia) - 1];
                }
              }
            }
            for (idx = 0; idx < nActiveConstr; idx++) {
              solution->lambda[idx] = -memspace->workspace_double[idx];
            }
          }
        }
        if ((solution->state != -7) ||
            (workingset->nActiveConstr > nVar_tmp_tmp)) {
          nActiveConstr = -1;
          minLambda = 0.0;
          n = (workingset->nWConstr[0] + workingset->nWConstr[1]) + 1;
          iyend = workingset->nActiveConstr;
          for (idx = n; idx <= iyend; idx++) {
            temp_tmp = solution->lambda[idx - 1];
            if (temp_tmp < minLambda) {
              minLambda = temp_tmp;
              nActiveConstr = idx - 1;
            }
          }
          if (nActiveConstr + 1 == 0) {
            solution->state = 1;
          } else {
            activeSetChangeID = -1;
            globalActiveConstrIdx = nActiveConstr + 1;
            subProblemChanged = true;
            iyend = workingset->Wid[nActiveConstr] - 1;
            workingset->isActiveConstr
                [(workingset->isActiveIdx[workingset->Wid[nActiveConstr] - 1] +
                  workingset->Wlocalidx[nActiveConstr]) -
                 2] = false;
            workingset->Wid[nActiveConstr] =
                workingset->Wid[workingset->nActiveConstr - 1];
            workingset->Wlocalidx[nActiveConstr] =
                workingset->Wlocalidx[workingset->nActiveConstr - 1];
            n = workingset->nVar;
            for (idx = 0; idx < n; idx++) {
              workingset->ATwset[idx + 5 * nActiveConstr] =
                  workingset->ATwset[idx + 5 * (workingset->nActiveConstr - 1)];
            }
            workingset->bwset[nActiveConstr] =
                workingset->bwset[workingset->nActiveConstr - 1];
            workingset->nActiveConstr--;
            workingset->nWConstr[iyend]--;
            solution->lambda[nActiveConstr] = 0.0;
          }
        } else {
          nActiveConstr = workingset->nActiveConstr;
          activeSetChangeID = 0;
          globalActiveConstrIdx = workingset->nActiveConstr;
          subProblemChanged = true;
          iyend = workingset->nActiveConstr - 1;
          idxRotGCol = workingset->Wid[iyend] - 1;
          workingset->isActiveConstr[(workingset->isActiveIdx[idxRotGCol] +
                                      workingset->Wlocalidx[iyend]) -
                                     2] = false;
          workingset->nActiveConstr--;
          workingset->nWConstr[idxRotGCol]--;
          solution->lambda[nActiveConstr - 1] = 0.0;
        }
        updateFval = false;
        stopping::checkStoppingAndUpdateFval(
            &activeSetChangeID, f, solution, memspace, objective, workingset,
            qrmanager, options_ObjectiveLimit, runTimeOptions_MaxIterations,
            updateFval);
      }
    } else {
      if (!updateFval) {
        solution->fstar = Objective::computeFval_ReuseHx(
            objective, memspace->workspace_double, f, solution->xstar);
      }
      exitg1 = 1;
    }
  } while (exitg1 == 0);
}

} // namespace qpactiveset
} // namespace coder
} // namespace optim
} // namespace coder

//
// File trailer for iterate.cpp
//
// [EOF]
//
