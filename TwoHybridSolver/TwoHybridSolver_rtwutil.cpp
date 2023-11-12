//
// File: TwoHybridSolver_rtwutil.cpp
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

// Include Files
#include "TwoHybridSolver_rtwutil.h"
#include "rt_nonfinite.h"
#include <string.h>

// Function Definitions
//
// Arguments    : int numerator
//                int denominator
// Return Type  : int
//
int div_nde_s32_floor(int numerator, int denominator)
{
  int b_numerator;
  if (((numerator < 0) != (denominator < 0)) &&
      (numerator % denominator != 0)) {
    b_numerator = -1;
  } else {
    b_numerator = 0;
  }
  return numerator / denominator + b_numerator;
}

//
// File trailer for TwoHybridSolver_rtwutil.cpp
//
// [EOF]
//
