//
// File: linearForm_.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

// Include Files
#include "linearForm_.h"
#include "rt_nonfinite.h"
#include <algorithm>
#include <cstring>
#include <string.h>

// Function Definitions
//
// Arguments    : boolean_T obj_hasLinear
//                int obj_nvar
//                double workspace[45]
//                const double H[16]
//                const double f[5]
//                const double x[5]
// Return Type  : void
//
namespace coder {
namespace optim {
namespace coder {
namespace qpactiveset {
namespace Objective {
void linearForm_(boolean_T obj_hasLinear, int obj_nvar, double workspace[45],
                 const double H[16], const double f[5], const double x[5])
{
  int ix;
  ix = 0;
  if (obj_hasLinear) {
    if (0 <= obj_nvar - 1) {
      std::copy(&f[0], &f[obj_nvar], &workspace[0]);
    }
    ix = 1;
  }
  if (obj_nvar != 0) {
    int i;
    if ((ix != 1) && (0 <= obj_nvar - 1)) {
      std::memset(&workspace[0], 0, obj_nvar * sizeof(double));
    }
    ix = 0;
    i = obj_nvar * (obj_nvar - 1) + 1;
    for (int iac{1}; obj_nvar < 0 ? iac >= i : iac <= i; iac += obj_nvar) {
      double c;
      int i1;
      c = 0.5 * x[ix];
      i1 = (iac + obj_nvar) - 1;
      for (int ia{iac}; ia <= i1; ia++) {
        int i2;
        i2 = ia - iac;
        workspace[i2] += H[ia - 1] * c;
      }
      ix++;
    }
  }
}

} // namespace Objective
} // namespace qpactiveset
} // namespace coder
} // namespace optim
} // namespace coder

//
// File trailer for linearForm_.cpp
//
// [EOF]
//
