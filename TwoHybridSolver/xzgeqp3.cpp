//
// File: xzgeqp3.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

// Include Files
#include "xzgeqp3.h"
#include "rt_nonfinite.h"
#include "xzlarf.h"
#include "xzlarfg.h"
#include <cstring>
#include <string.h>

// Function Definitions
//
// Arguments    : double A[81]
//                int m
//                int n
//                int nfxd
//                double tau[9]
// Return Type  : void
//
namespace coder {
namespace internal {
namespace reflapack {
void qrf(double A[81], int m, int n, int nfxd, double tau[9])
{
  double work[9];
  double atmp;
  std::memset(&work[0], 0, 9U * sizeof(double));
  for (int i{0}; i < nfxd; i++) {
    int ii;
    int mmi;
    ii = i * 9 + i;
    mmi = m - i;
    if (i + 1 < m) {
      atmp = A[ii];
      tau[i] = xzlarfg(mmi, &atmp, A, ii + 2);
      A[ii] = atmp;
    } else {
      tau[i] = 0.0;
    }
    if (i + 1 < n) {
      atmp = A[ii];
      A[ii] = 1.0;
      xzlarf(mmi, (n - i) - 1, ii + 1, tau[i], A, ii + 10, work);
      A[ii] = atmp;
    }
  }
}

} // namespace reflapack
} // namespace internal
} // namespace coder

//
// File trailer for xzgeqp3.cpp
//
// [EOF]
//
