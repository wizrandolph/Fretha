//
// File: _coder_TwoHybridSolver_mex.h
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

#ifndef _CODER_TWOHYBRIDSOLVER_MEX_H
#define _CODER_TWOHYBRIDSOLVER_MEX_H

// Include Files
#include "TwoHybridSolver/emlrt.h"
#include "TwoHybridSolver/mex.h"
#include "TwoHybridSolver/tmwtypes.h"

// Function Declarations
MEXFUNCTION_LINKAGE void mexFunction(int32_T nlhs, mxArray *plhs[],
                                     int32_T nrhs, const mxArray *prhs[]);

emlrtCTX mexFunctionCreateRootTLS();

void unsafe_TwoHybridSolver_mexFunction(int32_T nlhs, mxArray *plhs[2],
                                        int32_T nrhs, const mxArray *prhs[4]);

#endif
//
// File trailer for _coder_TwoHybridSolver_mex.h
//
// [EOF]
//
