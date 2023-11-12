//
// File: computeForwardDifferences.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

// Include Files
#include "computeForwardDifferences.h"
#include "TwoHybridSolver.h"
#include "TwoHybridSolver_internal_types.h"
#include "anonymous_function.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "rtGetInf.h"
#include <cmath>
#include <string.h>

// Function Definitions
//
// Arguments    : l_struct_T *obj
//                double fCurrent
//                double xk[4]
//                double gradf[5]
// Return Type  : boolean_T
//
namespace coder {
namespace optim {
namespace coder {
namespace utils {
namespace FiniteDifferences {
namespace internal {
boolean_T computeForwardDifferences(l_struct_T *obj, double fCurrent,
                                    double xk[4], double gradf[5])
{
  static double dv[4]{0.0, 0.0, 1.0, 1.0};
  int idx;
  boolean_T evalOK;
  boolean_T exitg1;
  dv[0U] = rtGetInf();
  dv[1U] = rtGetInf();
  evalOK = true;
  obj->numEvals = 0;
  idx = 0;
  exitg1 = false;
  while ((!exitg1) && (idx < 4)) {
    double d;
    double deltaX;
    double temp;
    double ubDiff;
    boolean_T guard1{false};
    boolean_T modifiedStep;
    deltaX = 1.4901161193847656E-8 *
             (1.0 - 2.0 * static_cast<double>(xk[idx] < 0.0)) *
             std::fmax(std::abs(xk[idx]), 1.0);
    if (obj->hasUB[idx]) {
      ubDiff = deltaX;
      modifiedStep = false;
      if ((xk[idx] >= 0.0) && (xk[idx] <= dv[idx]) &&
          ((xk[idx] + deltaX > dv[idx]) || (xk[idx] + deltaX < 0.0))) {
        ubDiff = -deltaX;
        modifiedStep = true;
        d = xk[idx] + -deltaX;
        if ((d > dv[idx]) || (d < 0.0)) {
          ubDiff = dv[idx] - xk[idx];
          if (xk[idx] <= ubDiff) {
            ubDiff = -xk[idx];
          }
        }
      }
      deltaX = ubDiff;
    } else if (obj->hasUB[idx]) {
      modifiedStep = false;
      if ((xk[idx] <= dv[idx]) && (xk[idx] + deltaX > dv[idx])) {
        deltaX = -deltaX;
        modifiedStep = true;
      }
    } else {
      modifiedStep = false;
      if ((xk[idx] >= 0.0) && (xk[idx] + deltaX < 0.0)) {
        deltaX = -deltaX;
        modifiedStep = true;
      }
    }
    temp = xk[idx];
    xk[idx] += deltaX;
    ubDiff =
        anon(obj->objfun.workspace.AAEST, obj->objfun.workspace.DDEST,
             obj->objfun.workspace.EACORRR, obj->objfun.workspace.EDCORRR, xk);
    evalOK = ((!std::isinf(ubDiff)) && (!std::isnan(ubDiff)));
    if (evalOK) {
      xk[idx] = temp;
    }
    obj->f_1 = ubDiff;
    obj->numEvals++;
    guard1 = false;
    if (!evalOK) {
      if (!modifiedStep) {
        deltaX = -deltaX;
        d = xk[idx] + deltaX;
        if ((d >= 0.0) && obj->hasUB[idx] && (d <= dv[idx])) {
          modifiedStep = true;
        } else {
          modifiedStep = false;
        }
        if ((!obj->hasBounds) || modifiedStep) {
          temp = xk[idx];
          xk[idx] = d;
          ubDiff = anon(
              obj->objfun.workspace.AAEST, obj->objfun.workspace.DDEST,
              obj->objfun.workspace.EACORRR, obj->objfun.workspace.EDCORRR, xk);
          evalOK = ((!std::isinf(ubDiff)) && (!std::isnan(ubDiff)));
          if (evalOK) {
            xk[idx] = temp;
          }
          obj->f_1 = ubDiff;
          obj->numEvals++;
        }
      }
      if (!evalOK) {
        exitg1 = true;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }
    if (guard1) {
      gradf[idx] = (obj->f_1 - fCurrent) / deltaX;
      idx++;
    }
  }
  return evalOK;
}

} // namespace internal
} // namespace FiniteDifferences
} // namespace utils
} // namespace coder
} // namespace optim
} // namespace coder

//
// File trailer for computeForwardDifferences.cpp
//
// [EOF]
//
