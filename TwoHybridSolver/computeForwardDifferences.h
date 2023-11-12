//
// File: computeForwardDifferences.h
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

#ifndef COMPUTEFORWARDDIFFERENCES_H
#define COMPUTEFORWARDDIFFERENCES_H

// Include Files
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

// Type Declarations
struct l_struct_T;

// Function Declarations
namespace coder {
namespace optim {
namespace coder {
namespace utils {
namespace FiniteDifferences {
namespace internal {
boolean_T computeForwardDifferences(l_struct_T *obj, double fCurrent,
                                    double xk[4], double gradf[5]);

}
} // namespace FiniteDifferences
} // namespace utils
} // namespace coder
} // namespace optim
} // namespace coder

#endif
//
// File trailer for computeForwardDifferences.h
//
// [EOF]
//
