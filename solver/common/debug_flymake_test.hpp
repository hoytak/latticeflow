#ifndef _DEBUG_FLYMAKE_TEST_H_
#define _DEBUG_FLYMAKE_TEST_H_

#ifdef EMACS_FLYMAKE

#ifndef _DEBUG_IN_TV_SOLVER_INCLUDES_

#include "../solvers.hpp"
#include "../kernels/kernels.hpp"
#include "../lattices/kernellattice.hpp"
#include "type_processing.hpp"

namespace latticeQBP {
  template RegPathPtr calculate2dTV<Star2d_4, long long>(size_t nx, size_t ny, 
                                                         double *function, 
                                                         const vector<double>& lambda);
};
#endif

#endif


#endif /* _DEBUG_FLYMAKE_TEST_H_ */
