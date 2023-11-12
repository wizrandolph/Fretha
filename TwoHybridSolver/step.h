//
// File: step.h
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

#ifndef STEP_H
#define STEP_H

// Include Files
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

// Type Declarations
struct i_struct_T;

struct b_struct_T;

struct d_struct_T;

struct j_struct_T;

struct e_struct_T;

struct f_struct_T;

struct struct_T;

struct g_struct_T;

// Function Declarations
namespace coder {
namespace optim {
namespace coder {
namespace fminconsqp {
boolean_T step(int *STEP_TYPE, double Hessian[16], i_struct_T *TrialState,
               b_struct_T *MeritFunction, d_struct_T *memspace,
               j_struct_T *WorkingSet, e_struct_T *b_QRManager,
               f_struct_T *b_CholManager, struct_T *QPObjective,
               g_struct_T *qpoptions);

}
} // namespace coder
} // namespace optim
} // namespace coder

#endif
//
// File trailer for step.h
//
// [EOF]
//
