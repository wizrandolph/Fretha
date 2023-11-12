//
// File: BFGSUpdate.h
//
// MATLAB Coder version            : 5.2
// C/C++ source code generated on  : 16-Mar-2023 17:31:52
//

#ifndef BFGSUPDATE_H
#define BFGSUPDATE_H

// Include Files
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
namespace optim {
namespace coder {
namespace fminconsqp {
boolean_T BFGSUpdate(int nvar, double Bk[16], const double sk[5], double yk[5],
                     double workspace[45]);

}
} // namespace coder
} // namespace optim
} // namespace coder

#endif
//
// File trailer for BFGSUpdate.h
//
// [EOF]
//
