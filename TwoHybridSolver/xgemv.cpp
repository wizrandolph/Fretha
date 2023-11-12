//
// File: xgemv.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

// Include Files
#include "xgemv.h"
#include "TwoHybridSolver_rtwutil.h"
#include "rt_nonfinite.h"
#include <cstring>
#include <string.h>

// Function Definitions
//
// Arguments    : int m
//                int n
//                const double A[81]
//                const double x[5]
//                double y[45]
// Return Type  : void
//
namespace coder {
namespace internal {
namespace blas {
void xgemv(int m, int n, const double A[81], const double x[5], double y[45])
{
  if (m != 0) {
    int i;
    if (0 <= n - 1) {
      std::memset(&y[0], 0, n * sizeof(double));
    }
    i = 9 * (n - 1) + 1;
    for (int iac{1}; iac <= i; iac += 9) {
      double c;
      int i1;
      c = 0.0;
      i1 = (iac + m) - 1;
      for (int ia{iac}; ia <= i1; ia++) {
        c += A[ia - 1] * x[ia - iac];
      }
      i1 = div_nde_s32_floor(iac - 1, 9);
      y[i1] += c;
    }
  }
}

} // namespace blas
} // namespace internal
} // namespace coder

//
// File trailer for xgemv.cpp
//
// [EOF]
//
