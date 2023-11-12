//
// File: computeFval_ReuseHx.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

// Include Files
#include "computeFval_ReuseHx.h"
#include "TwoHybridSolver_internal_types.h"
#include "rt_nonfinite.h"
#include <algorithm>
#include <string.h>

// Function Definitions
//
// Arguments    : const struct_T *obj
//                double workspace[45]
//                const double f[5]
//                const double x[5]
// Return Type  : double
//
namespace coder {
namespace optim {
namespace coder {
namespace qpactiveset {
namespace Objective {
double computeFval_ReuseHx(const struct_T *obj, double workspace[45],
                           const double f[5], const double x[5])
{
  double val;
  switch (obj->objtype) {
  case 5:
    val = obj->gammaScalar * x[obj->nvar - 1];
    break;
  case 3: {
    if (obj->hasLinear) {
      int idx;
      int ixlast;
      ixlast = obj->nvar;
      for (idx = 0; idx < ixlast; idx++) {
        workspace[idx] = 0.5 * obj->Hx[idx] + f[idx];
      }
      val = 0.0;
      if (obj->nvar >= 1) {
        ixlast = obj->nvar;
        for (idx = 0; idx < ixlast; idx++) {
          val += x[idx] * workspace[idx];
        }
      }
    } else {
      val = 0.0;
      if (obj->nvar >= 1) {
        int ixlast;
        ixlast = obj->nvar;
        for (int idx{0}; idx < ixlast; idx++) {
          val += x[idx] * obj->Hx[idx];
        }
      }
      val *= 0.5;
    }
  } break;
  default: {
    if (obj->hasLinear) {
      int ixlast;
      ixlast = obj->nvar;
      if (0 <= ixlast - 1) {
        std::copy(&f[0], &f[ixlast], &workspace[0]);
      }
      ixlast = 3 - obj->nvar;
      for (int idx{0}; idx <= ixlast; idx++) {
        workspace[obj->nvar + idx] = obj->rho;
      }
      workspace[0] += 0.5 * obj->Hx[0];
      workspace[1] += 0.5 * obj->Hx[1];
      workspace[2] += 0.5 * obj->Hx[2];
      workspace[3] += 0.5 * obj->Hx[3];
      val =
          ((x[0] * workspace[0] + x[1] * workspace[1]) + x[2] * workspace[2]) +
          x[3] * workspace[3];
    } else {
      int ixlast;
      val =
          0.5 * (((x[0] * obj->Hx[0] + x[1] * obj->Hx[1]) + x[2] * obj->Hx[2]) +
                 x[3] * obj->Hx[3]);
      ixlast = obj->nvar + 1;
      for (int idx{ixlast}; idx < 5; idx++) {
        val += x[idx - 1] * obj->rho;
      }
    }
  } break;
  }
  return val;
}

} // namespace Objective
} // namespace qpactiveset
} // namespace coder
} // namespace optim
} // namespace coder

//
// File trailer for computeFval_ReuseHx.cpp
//
// [EOF]
//
