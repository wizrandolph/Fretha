//
// File: driver.h
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

#ifndef DRIVER_H
#define DRIVER_H

// Include Files
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

// Type Declarations
struct i_struct_T;

struct b_struct_T;

struct k_struct_T;

struct l_struct_T;

struct d_struct_T;

struct j_struct_T;

struct e_struct_T;

struct struct_T;

struct f_struct_T;

// Function Declarations
namespace coder {
namespace optim {
namespace coder {
namespace fminconsqp {
void driver(i_struct_T *TrialState, b_struct_T *MeritFunction,
            const k_struct_T *FcnEvaluator, l_struct_T *FiniteDifferences,
            d_struct_T *memspace, j_struct_T *WorkingSet,
            e_struct_T *b_QRManager, struct_T *QPObjective, double Hessian[16],
            f_struct_T *b_CholManager);

}
} // namespace coder
} // namespace optim
} // namespace coder

#endif
//
// File trailer for driver.h
//
// [EOF]
//
