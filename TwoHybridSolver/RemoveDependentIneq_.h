//
// File: RemoveDependentIneq_.h
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

#ifndef REMOVEDEPENDENTINEQ__H
#define REMOVEDEPENDENTINEQ__H

// Include Files
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

// Type Declarations
struct j_struct_T;

struct e_struct_T;

struct d_struct_T;

// Function Declarations
namespace coder {
namespace optim {
namespace coder {
namespace qpactiveset {
namespace initialize {
void RemoveDependentIneq_(j_struct_T *workingset, e_struct_T *qrmanager,
                          d_struct_T *memspace, double tolfactor);

}
} // namespace qpactiveset
} // namespace coder
} // namespace optim
} // namespace coder

#endif
//
// File trailer for RemoveDependentIneq_.h
//
// [EOF]
//
