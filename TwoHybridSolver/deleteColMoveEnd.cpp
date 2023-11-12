//
// File: deleteColMoveEnd.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

// Include Files
#include "deleteColMoveEnd.h"
#include "TwoHybridSolver_internal_types.h"
#include "rt_nonfinite.h"
#include "xrotg.h"
#include <string.h>

// Function Definitions
//
// Arguments    : e_struct_T *obj
//                int idx
// Return Type  : void
//
namespace coder {
namespace optim {
namespace coder {
namespace QRManager {
void deleteColMoveEnd(e_struct_T *obj, int idx)
{
  double c;
  double s;
  double temp_tmp;
  int i;
  if (obj->usedPivoting) {
    i = 1;
    while ((i <= obj->ncols) && (obj->jpvt[i - 1] != idx)) {
      i++;
    }
    idx = i;
  }
  if (idx >= obj->ncols) {
    obj->ncols--;
  } else {
    int b_i;
    int k;
    int u0;
    obj->jpvt[idx - 1] = obj->jpvt[obj->ncols - 1];
    b_i = obj->minRowCol;
    for (k = 0; k < b_i; k++) {
      obj->QR[k + 9 * (idx - 1)] = obj->QR[k + 9 * (obj->ncols - 1)];
    }
    obj->ncols--;
    u0 = obj->mrows;
    i = obj->ncols;
    if (u0 < i) {
      i = u0;
    }
    obj->minRowCol = i;
    if (idx < obj->mrows) {
      double c_temp_tmp;
      int QRk0;
      int b_k;
      int b_temp_tmp;
      int endIdx;
      int n;
      u0 = obj->mrows - 1;
      endIdx = obj->ncols;
      if (u0 < endIdx) {
        endIdx = u0;
      }
      k = endIdx;
      i = 9 * (idx - 1);
      while (k >= idx) {
        b_i = k + i;
        temp_tmp = obj->QR[b_i];
        internal::blas::xrotg(&obj->QR[(k + i) - 1], &temp_tmp, &c, &s);
        obj->QR[b_i] = temp_tmp;
        b_i = 9 * (k - 1);
        obj->QR[k + b_i] = 0.0;
        QRk0 = k + 9 * idx;
        n = obj->ncols - idx;
        if (n >= 1) {
          for (b_k = 0; b_k < n; b_k++) {
            b_temp_tmp = QRk0 + b_k * 9;
            temp_tmp = obj->QR[b_temp_tmp];
            c_temp_tmp = obj->QR[b_temp_tmp - 1];
            obj->QR[b_temp_tmp] = c * temp_tmp - s * c_temp_tmp;
            obj->QR[b_temp_tmp - 1] = c * c_temp_tmp + s * temp_tmp;
          }
        }
        n = obj->mrows;
        if (obj->mrows >= 1) {
          for (b_k = 0; b_k < n; b_k++) {
            b_temp_tmp = b_i + b_k;
            temp_tmp = obj->Q[b_temp_tmp + 9];
            c_temp_tmp = obj->Q[b_temp_tmp];
            obj->Q[b_temp_tmp + 9] = c * temp_tmp - s * c_temp_tmp;
            obj->Q[b_temp_tmp] = c * c_temp_tmp + s * temp_tmp;
          }
        }
        k--;
      }
      b_i = idx + 1;
      for (k = b_i; k <= endIdx; k++) {
        u0 = 9 * (k - 1);
        i = k + u0;
        temp_tmp = obj->QR[i];
        internal::blas::xrotg(&obj->QR[(k + u0) - 1], &temp_tmp, &c, &s);
        obj->QR[i] = temp_tmp;
        QRk0 = k * 10;
        n = obj->ncols - k;
        if (n >= 1) {
          for (b_k = 0; b_k < n; b_k++) {
            b_temp_tmp = QRk0 + b_k * 9;
            temp_tmp = obj->QR[b_temp_tmp];
            c_temp_tmp = obj->QR[b_temp_tmp - 1];
            obj->QR[b_temp_tmp] = c * temp_tmp - s * c_temp_tmp;
            obj->QR[b_temp_tmp - 1] = c * c_temp_tmp + s * temp_tmp;
          }
        }
        n = obj->mrows;
        if (obj->mrows >= 1) {
          for (b_k = 0; b_k < n; b_k++) {
            b_temp_tmp = u0 + b_k;
            temp_tmp = obj->Q[b_temp_tmp + 9];
            c_temp_tmp = obj->Q[b_temp_tmp];
            obj->Q[b_temp_tmp + 9] = c * temp_tmp - s * c_temp_tmp;
            obj->Q[b_temp_tmp] = c * c_temp_tmp + s * temp_tmp;
          }
        }
      }
    }
  }
}

} // namespace QRManager
} // namespace coder
} // namespace optim
} // namespace coder

//
// File trailer for deleteColMoveEnd.cpp
//
// [EOF]
//
