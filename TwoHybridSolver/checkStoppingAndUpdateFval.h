//
// File: checkStoppingAndUpdateFval.h
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

#ifndef CHECKSTOPPINGANDUPDATEFVAL_H
#define CHECKSTOPPINGANDUPDATEFVAL_H

// Include Files
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

// Type Declarations
struct i_struct_T;

struct d_struct_T;

struct struct_T;

struct j_struct_T;

struct e_struct_T;

// Function Declarations
namespace coder {
namespace optim {
namespace coder {
namespace qpactiveset {
namespace stopping {
void checkStoppingAndUpdateFval(int *activeSetChangeID, const double f[5],
                                i_struct_T *solution, d_struct_T *memspace,
                                const struct_T *objective,
                                const j_struct_T *workingset,
                                e_struct_T *qrmanager,
                                double options_ObjectiveLimit,
                                int runTimeOptions_MaxIterations,
                                boolean_T updateFval);

}
} // namespace qpactiveset
} // namespace coder
} // namespace optim
} // namespace coder

#endif
//
// File trailer for checkStoppingAndUpdateFval.h
//
// [EOF]
//
