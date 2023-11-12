//
// File: setProblemType.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

// Include Files
#include "setProblemType.h"
#include "TwoHybridSolver_internal_types.h"
#include "rt_nonfinite.h"
#include <cstring>
#include <string.h>

// Function Definitions
//
// Arguments    : j_struct_T *obj
//                int PROBLEM_TYPE
// Return Type  : void
//
namespace coder {
namespace optim {
namespace coder {
namespace qpactiveset {
namespace WorkingSet {
void setProblemType(j_struct_T *obj, int PROBLEM_TYPE)
{
  switch (PROBLEM_TYPE) {
  case 3: {
    int i;
    obj->nVar = 4;
    obj->mConstr = 6;
    if (obj->nWConstr[4] > 0) {
      obj->isActiveConstr[4] = obj->isActiveConstr[obj->isActiveIdx[4] - 1];
      obj->isActiveConstr[5] = obj->isActiveConstr[obj->isActiveIdx[4]];
    }
    for (i = 0; i < 5; i++) {
      obj->sizes[i] = obj->sizesNormal[i];
    }
    for (i = 0; i < 6; i++) {
      obj->isActiveIdx[i] = obj->isActiveIdxNormal[i];
    }
  } break;
  case 1: {
    int i;
    obj->nVar = 5;
    obj->mConstr = 7;
    for (i = 0; i < 5; i++) {
      obj->sizes[i] = obj->sizesPhaseOne[i];
    }
    for (i = 0; i < 6; i++) {
      obj->isActiveIdx[i] = obj->isActiveIdxPhaseOne[i];
    }
    obj->indexLB[4] = 5;
    obj->lb[4] = 1.0E-5;
    i = obj->nActiveConstr;
    for (int colOffsetATw{1}; colOffsetATw <= i; colOffsetATw++) {
      obj->ATwset[5 * (colOffsetATw - 1) + 4] = -1.0;
    }
    if (obj->nWConstr[4] > 0) {
      obj->isActiveConstr[6] = obj->isActiveConstr[5];
      obj->isActiveConstr[7] = obj->isActiveConstr[6];
    }
    obj->isActiveConstr[4] = false;
  } break;
  case 2: {
    int i;
    obj->nVar = 4;
    obj->mConstr = 8;
    for (i = 0; i < 5; i++) {
      obj->sizes[i] = obj->sizesRegularized[i];
    }
    if (obj->probType != 4) {
      if (obj->nWConstr[4] > 0) {
        obj->isActiveConstr[5] = obj->isActiveConstr[obj->isActiveIdx[4] - 1];
        obj->isActiveConstr[6] = obj->isActiveConstr[obj->isActiveIdx[4]];
      }
      i = obj->nActiveConstr;
      for (int idx_col{1}; idx_col <= i; idx_col++) {
        int colOffsetATw;
        colOffsetATw = 5 * (idx_col - 1) - 1;
        switch (obj->Wid[idx_col - 1]) {
        case 3: {
          int i1;
          i1 = obj->Wlocalidx[idx_col - 1] + 3;
          if (5 <= i1) {
            std::memset(&obj->ATwset[colOffsetATw + 5], 0,
                        (((i1 + colOffsetATw) - colOffsetATw) + -4) *
                            sizeof(double));
          }
          obj->ATwset[(obj->Wlocalidx[idx_col - 1] + colOffsetATw) + 4] = -1.0;
          i1 = obj->Wlocalidx[idx_col - 1] + 5;
          if (i1 <= 4) {
            std::memset(&obj->ATwset[i1 + colOffsetATw], 0,
                        (((colOffsetATw - i1) - colOffsetATw) + 5) *
                            sizeof(double));
          }
        } break;
        }
      }
    }
    for (i = 0; i < 6; i++) {
      obj->isActiveIdx[i] = obj->isActiveIdxRegularized[i];
    }
  } break;
  default: {
    int i;
    obj->nVar = 5;
    obj->mConstr = 9;
    for (i = 0; i < 5; i++) {
      obj->sizes[i] = obj->sizesRegPhaseOne[i];
    }
    for (i = 0; i < 6; i++) {
      obj->isActiveIdx[i] = obj->isActiveIdxRegPhaseOne[i];
    }
    obj->indexLB[4] = 5;
    obj->lb[4] = 1.0E-5;
    i = obj->nActiveConstr;
    for (int colOffsetATw{1}; colOffsetATw <= i; colOffsetATw++) {
      obj->ATwset[5 * (colOffsetATw - 1) + 4] = -1.0;
    }
    if (obj->nWConstr[4] > 0) {
      obj->isActiveConstr[6] = obj->isActiveConstr[5];
      obj->isActiveConstr[7] = obj->isActiveConstr[6];
    }
    obj->isActiveConstr[4] = false;
  } break;
  }
  obj->probType = PROBLEM_TYPE;
}

} // namespace WorkingSet
} // namespace qpactiveset
} // namespace coder
} // namespace optim
} // namespace coder

//
// File trailer for setProblemType.cpp
//
// [EOF]
//
