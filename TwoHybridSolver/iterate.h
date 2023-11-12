//
// File: iterate.h
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

#ifndef ITERATE_H
#define ITERATE_H

// Include Files
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

// Type Declarations
struct i_struct_T;

struct d_struct_T;

struct j_struct_T;

struct e_struct_T;

struct f_struct_T;

struct struct_T;

// Function Declarations
namespace coder {
namespace optim {
namespace coder {
namespace qpactiveset {
void iterate(const double H[16], const double f[5], i_struct_T *solution,
             d_struct_T *memspace, j_struct_T *workingset,
             e_struct_T *qrmanager, f_struct_T *cholmanager,
             struct_T *objective, const char options_SolverName[7],
             double options_StepTolerance, double options_ObjectiveLimit,
             int runTimeOptions_MaxIterations);

}
} // namespace coder
} // namespace optim
} // namespace coder

#endif
//
// File trailer for iterate.h
//
// [EOF]
//
