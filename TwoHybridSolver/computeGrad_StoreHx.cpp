//
// File: computeGrad_StoreHx.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

// Include Files
#include "computeGrad_StoreHx.h"
#include "TwoHybridSolver_internal_types.h"
#include "rt_nonfinite.h"
#include <algorithm>
#include <cstring>
#include <string.h>

// Function Definitions
//
// Arguments    : struct_T *obj
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
void computeGrad_StoreHx(struct_T *obj, const double H[16], const double f[5],
                         const double x[5])
{
  switch (obj->objtype) {
  case 5: {
    int i;
    i = obj->nvar;
    if (0 <= i - 2) {
      std::memset(&obj->grad[0], 0, (i + -1) * sizeof(double));
    }
    obj->grad[obj->nvar - 1] = obj->gammaScalar;
  } break;
  case 3: {
    int i;
    int iy;
    int lda;
    iy = obj->nvar - 1;
    lda = obj->nvar;
    if (obj->nvar != 0) {
      int ix;
      if (0 <= iy) {
        std::memset(&obj->Hx[0], 0, (iy + 1) * sizeof(double));
      }
      ix = 0;
      i = obj->nvar * (obj->nvar - 1) + 1;
      for (int iac{1}; lda < 0 ? iac >= i : iac <= i; iac += lda) {
        int i1;
        i1 = iac + iy;
        for (int ia{iac}; ia <= i1; ia++) {
          int i2;
          i2 = ia - iac;
          obj->Hx[i2] += H[ia - 1] * x[ix];
        }
        ix++;
      }
    }
    i = obj->nvar;
    if (0 <= i - 1) {
      std::copy(&obj->Hx[0], &obj->Hx[i], &obj->grad[0]);
    }
    if (obj->hasLinear && (obj->nvar >= 1)) {
      i = obj->nvar - 1;
      for (lda = 0; lda <= i; lda++) {
        obj->grad[lda] += f[lda];
      }
    }
  } break;
  default: {
    int i;
    int i1;
    int iy;
    int lda;
    iy = obj->nvar - 1;
    lda = obj->nvar;
    if (obj->nvar != 0) {
      int ix;
      if (0 <= iy) {
        std::memset(&obj->Hx[0], 0, (iy + 1) * sizeof(double));
      }
      ix = 0;
      i = obj->nvar * (obj->nvar - 1) + 1;
      for (int iac{1}; lda < 0 ? iac >= i : iac <= i; iac += lda) {
        i1 = iac + iy;
        for (int ia{iac}; ia <= i1; ia++) {
          int i2;
          i2 = ia - iac;
          obj->Hx[i2] += H[ia - 1] * x[ix];
        }
        ix++;
      }
    }
    i = obj->nvar + 1;
    for (iy = i; iy < 5; iy++) {
      obj->Hx[iy - 1] = obj->beta * x[iy - 1];
    }
    obj->grad[0] = obj->Hx[0];
    obj->grad[1] = obj->Hx[1];
    obj->grad[2] = obj->Hx[2];
    obj->grad[3] = obj->Hx[3];
    if (obj->hasLinear && (obj->nvar >= 1)) {
      i = obj->nvar - 1;
      for (lda = 0; lda <= i; lda++) {
        obj->grad[lda] += f[lda];
      }
    }
    if (4 - obj->nvar >= 1) {
      iy = obj->nvar;
      i = 3 - obj->nvar;
      for (lda = 0; lda <= i; lda++) {
        i1 = iy + lda;
        obj->grad[i1] += obj->rho;
      }
    }
  } break;
  }
}

} // namespace Objective
} // namespace qpactiveset
} // namespace coder
} // namespace optim
} // namespace coder

//
// File trailer for computeGrad_StoreHx.cpp
//
// [EOF]
//
