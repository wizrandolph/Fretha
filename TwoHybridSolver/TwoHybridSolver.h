//
// File: TwoHybridSolver.h
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

#ifndef TWOHYBRIDSOLVER_H
#define TWOHYBRIDSOLVER_H

// Include Files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
extern void TwoHybridSolver(const coder::array<double, 2U> &AAEST,
                            const coder::array<double, 2U> &DDEST,
                            const coder::array<double, 2U> &EACORRR,
                            const coder::array<double, 2U> &EDCORRR,
                            double x[4], double *y);

double anon(const coder::array<double, 2U> &AAEST,
            const coder::array<double, 2U> &DDEST,
            const coder::array<double, 2U> &EACORRR,
            const coder::array<double, 2U> &EDCORRR, const double X[4]);

#endif
//
// File trailer for TwoHybridSolver.h
//
// [EOF]
//
