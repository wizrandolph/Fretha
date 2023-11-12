//
// File: solve.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

// Include Files
#include "solve.h"
#include "TwoHybridSolver_internal_types.h"
#include "rt_nonfinite.h"
#include <string.h>

// Function Definitions
//
// Arguments    : const f_struct_T *obj
//                double rhs[5]
// Return Type  : void
//
namespace coder {
namespace optim {
namespace coder {
namespace CholManager {
void solve(const f_struct_T *obj, double rhs[5])
{
  int n_tmp;
  n_tmp = obj->ndims;
  if (obj->ndims != 0) {
    int i;
    int j;
    int jA;
    for (j = 0; j < n_tmp; j++) {
      double temp;
      jA = j * 9;
      temp = rhs[j];
      for (i = 0; i < j; i++) {
        temp -= obj->FMat[jA + i] * rhs[i];
      }
      rhs[j] = temp / obj->FMat[jA + j];
    }
    for (j = n_tmp; j >= 1; j--) {
      jA = (j + (j - 1) * 9) - 1;
      rhs[j - 1] /= obj->FMat[jA];
      for (i = 0; i <= j - 2; i++) {
        int ix;
        ix = (j - i) - 2;
        rhs[ix] -= rhs[j - 1] * obj->FMat[(jA - i) - 1];
      }
    }
  }
}

} // namespace CholManager
} // namespace coder
} // namespace optim
} // namespace coder

//
// File trailer for solve.cpp
//
// [EOF]
//
