//
// File: test_exit.h
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

#ifndef TEST_EXIT_H
#define TEST_EXIT_H

// Include Files
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

// Type Declarations
struct b_struct_T;

struct j_struct_T;

struct i_struct_T;

struct c_struct_T;

struct d_struct_T;

struct e_struct_T;

// Function Declarations
namespace coder {
namespace optim {
namespace coder {
namespace fminconsqp {
void b_test_exit(c_struct_T *Flags, d_struct_T *memspace,
                 b_struct_T *MeritFunction, const j_struct_T *WorkingSet,
                 i_struct_T *TrialState, e_struct_T *b_QRManager);

void test_exit(b_struct_T *MeritFunction, const j_struct_T *WorkingSet,
               i_struct_T *TrialState, boolean_T *Flags_gradOK,
               boolean_T *Flags_fevalOK, boolean_T *Flags_done,
               boolean_T *Flags_stepAccepted, boolean_T *Flags_failedLineSearch,
               int *Flags_stepType);

} // namespace fminconsqp
} // namespace coder
} // namespace optim
} // namespace coder

#endif
//
// File trailer for test_exit.h
//
// [EOF]
//
