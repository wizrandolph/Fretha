//
// File: fullColLDL2_.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

// Include Files
#include "fullColLDL2_.h"
#include "TwoHybridSolver_internal_types.h"
#include "rt_nonfinite.h"
#include <string.h>

// Function Definitions
//
// Arguments    : f_struct_T *obj
//                int NColsRemain
// Return Type  : void
//
namespace coder {
namespace optim {
namespace coder {
namespace DynamicRegCholManager {
void fullColLDL2_(f_struct_T *obj, int NColsRemain)
{
  for (int k{0}; k < NColsRemain; k++) {
    double alpha1;
    double y;
    int LD_diagOffset;
    int i;
    int jA;
    int subMatrixDim;
    LD_diagOffset = 10 * k;
    alpha1 = -1.0 / obj->FMat[LD_diagOffset];
    subMatrixDim = (NColsRemain - k) - 2;
    for (jA = 0; jA <= subMatrixDim; jA++) {
      obj->workspace_ = obj->FMat[(LD_diagOffset + jA) + 1];
    }
    y = obj->workspace_;
    if (!(alpha1 == 0.0)) {
      jA = LD_diagOffset;
      for (int j{0}; j <= subMatrixDim; j++) {
        if (y != 0.0) {
          double temp;
          int i1;
          temp = y * alpha1;
          i = jA + 11;
          i1 = subMatrixDim + jA;
          for (int ijA{i}; ijA <= i1 + 11; ijA++) {
            obj->FMat[ijA - 1] += obj->workspace_ * temp;
          }
        }
        jA += 9;
      }
    }
    for (jA = 0; jA <= subMatrixDim; jA++) {
      i = (LD_diagOffset + jA) + 1;
      obj->FMat[i] /= obj->FMat[LD_diagOffset];
    }
  }
}

} // namespace DynamicRegCholManager
} // namespace coder
} // namespace optim
} // namespace coder

//
// File trailer for fullColLDL2_.cpp
//
// [EOF]
//
