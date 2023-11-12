//
// File: sortLambdaQP.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

// Include Files
#include "sortLambdaQP.h"
#include "rt_nonfinite.h"
#include <string.h>

// Function Definitions
//
// Arguments    : double lambda[9]
//                int WorkingSet_nActiveConstr
//                const int WorkingSet_sizes[5]
//                const int WorkingSet_isActiveIdx[6]
//                const int WorkingSet_Wid[9]
//                const int WorkingSet_Wlocalidx[9]
//                double workspace[45]
// Return Type  : void
//
namespace coder {
namespace optim {
namespace coder {
namespace qpactiveset {
namespace parseoutput {
void sortLambdaQP(double lambda[9], int WorkingSet_nActiveConstr,
                  const int WorkingSet_sizes[5],
                  const int WorkingSet_isActiveIdx[6],
                  const int WorkingSet_Wid[9],
                  const int WorkingSet_Wlocalidx[9], double workspace[45])
{
  if (WorkingSet_nActiveConstr != 0) {
    int idx;
    int mAll;
    mAll = WorkingSet_sizes[3] + 1;
    for (idx = 0; idx <= mAll; idx++) {
      workspace[idx] = lambda[idx];
      lambda[idx] = 0.0;
    }
    mAll = 0;
    idx = 0;
    while ((idx + 1 <= WorkingSet_nActiveConstr) &&
           (WorkingSet_Wid[idx] <= 2)) {
      lambda[WorkingSet_Wlocalidx[idx] - 1] = workspace[mAll];
      mAll++;
      idx++;
    }
    while (idx + 1 <= WorkingSet_nActiveConstr) {
      int idxOffset;
      switch (WorkingSet_Wid[idx]) {
      case 3:
        idxOffset = 1;
        break;
      case 4:
        idxOffset = 1;
        break;
      default:
        idxOffset = WorkingSet_isActiveIdx[4];
        break;
      }
      lambda[(idxOffset + WorkingSet_Wlocalidx[idx]) - 2] = workspace[mAll];
      mAll++;
      idx++;
    }
  }
}

} // namespace parseoutput
} // namespace qpactiveset
} // namespace coder
} // namespace optim
} // namespace coder

//
// File trailer for sortLambdaQP.cpp
//
// [EOF]
//
