//
// File: computeFval.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

// Include Files
#include "computeFval.h"
#include "TwoHybridSolver_internal_types.h"
#include "linearForm_.h"
#include "rt_nonfinite.h"
#include <string.h>

// Function Definitions
//
// Arguments    : const struct_T *obj
//                double workspace[45]
//                const double H[16]
//                const double f[5]
//                const double x[5]
// Return Type  : double
//
namespace coder {
namespace optim {
namespace coder {
namespace qpactiveset {
namespace Objective {
double computeFval(const struct_T *obj, double workspace[45],
                   const double H[16], const double f[5], const double x[5])
{
  double val;
  switch (obj->objtype) {
  case 5:
    val = obj->gammaScalar * x[obj->nvar - 1];
    break;
  case 3: {
    linearForm_(obj->hasLinear, obj->nvar, workspace, H, f, x);
    val = 0.0;
    if (obj->nvar >= 1) {
      int ixlast;
      ixlast = obj->nvar;
      for (int idx{0}; idx < ixlast; idx++) {
        val += x[idx] * workspace[idx];
      }
    }
  } break;
  default: {
    int ixlast;
    linearForm_(obj->hasLinear, obj->nvar, workspace, H, f, x);
    ixlast = obj->nvar + 1;
    for (int idx{ixlast}; idx < 5; idx++) {
      workspace[idx - 1] = 0.5 * obj->beta * x[idx - 1] + obj->rho;
    }
    val = ((x[0] * workspace[0] + x[1] * workspace[1]) + x[2] * workspace[2]) +
          x[3] * workspace[3];
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
// File trailer for computeFval.cpp
//
// [EOF]
//
