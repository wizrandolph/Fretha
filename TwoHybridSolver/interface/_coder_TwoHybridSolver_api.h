//
// File: _coder_TwoHybridSolver_api.h
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

#ifndef _CODER_TWOHYBRIDSOLVER_API_H
#define _CODER_TWOHYBRIDSOLVER_API_H

// Include Files
#include "coder_array_mex.h"
#include "../emlrt.h"
#include "../tmwtypes.h"
#include <algorithm>
#include <cstring>

// Variable Declarations
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

// Function Declarations
void TwoHybridSolver(coder::array<real_T, 2U> *AAEST,
                     coder::array<real_T, 2U> *DDEST,
                     coder::array<real_T, 2U> *EACORRR,
                     coder::array<real_T, 2U> *EDCORRR, real_T x[4], real_T *y);

void TwoHybridSolver_api(const mxArray *const prhs[4], int32_T nlhs,
                         const mxArray *plhs[2]);

void TwoHybridSolver_atexit();

void TwoHybridSolver_initialize();

void TwoHybridSolver_terminate();

void TwoHybridSolver_xil_shutdown();

void TwoHybridSolver_xil_terminate();

#endif
//
// File trailer for _coder_TwoHybridSolver_api.h
//
// [EOF]
//
