//
// File: xgemm.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

// Include Files
#include "xgemm.h"
#include "rt_nonfinite.h"
#include <cstring>
#include <string.h>

// Function Definitions
//
// Arguments    : int m
//                int n
//                int k
//                const double A[16]
//                int lda
//                const double B[81]
//                int ib0
//                double C[45]
// Return Type  : void
//
namespace coder {
namespace internal {
namespace blas {
void xgemm(int m, int n, int k, const double A[16], int lda, const double B[81],
           int ib0, double C[45])
{
  if ((m != 0) && (n != 0)) {
    int br;
    int cr;
    int i;
    int i1;
    int lastColC;
    br = ib0;
    lastColC = 9 * (n - 1);
    for (cr = 0; cr <= lastColC; cr += 9) {
      i = cr + 1;
      i1 = cr + m;
      if (i <= i1) {
        std::memset(&C[i + -1], 0, ((i1 - i) + 1) * sizeof(double));
      }
    }
    for (cr = 0; cr <= lastColC; cr += 9) {
      int ar;
      ar = -1;
      i = br + k;
      for (int ib{br}; ib < i; ib++) {
        int i2;
        i1 = cr + 1;
        i2 = cr + m;
        for (int ic{i1}; ic <= i2; ic++) {
          C[ic - 1] += B[ib - 1] * A[(ar + ic) - cr];
        }
        ar += lda;
      }
      br += 9;
    }
  }
}

//
// Arguments    : int m
//                int n
//                int k
//                const double A[81]
//                int ia0
//                const double B[45]
//                double C[81]
// Return Type  : void
//
void xgemm(int m, int n, int k, const double A[81], int ia0, const double B[45],
           double C[81])
{
  if ((m != 0) && (n != 0)) {
    int br;
    int cr;
    int i;
    int i1;
    int lastColC;
    lastColC = 9 * (n - 1);
    for (cr = 0; cr <= lastColC; cr += 9) {
      i = cr + 1;
      i1 = cr + m;
      if (i <= i1) {
        std::memset(&C[i + -1], 0, ((i1 - i) + 1) * sizeof(double));
      }
    }
    br = -1;
    for (cr = 0; cr <= lastColC; cr += 9) {
      int ar;
      ar = ia0;
      i = cr + 1;
      i1 = cr + m;
      for (int ic{i}; ic <= i1; ic++) {
        double temp;
        temp = 0.0;
        for (int w{0}; w < k; w++) {
          temp += A[(w + ar) - 1] * B[(w + br) + 1];
        }
        C[ic - 1] += temp;
        ar += 9;
      }
      br += 9;
    }
  }
}

} // namespace blas
} // namespace internal
} // namespace coder

//
// File trailer for xgemm.cpp
//
// [EOF]
//
