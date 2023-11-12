//
// File: xgemm.h
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

#ifndef XGEMM_H
#define XGEMM_H

// Include Files
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
namespace internal {
namespace blas {
void xgemm(int m, int n, int k, const double A[16], int lda, const double B[81],
           int ib0, double C[45]);

void xgemm(int m, int n, int k, const double A[81], int ia0, const double B[45],
           double C[81]);

} // namespace blas
} // namespace internal
} // namespace coder

#endif
//
// File trailer for xgemm.h
//
// [EOF]
//
