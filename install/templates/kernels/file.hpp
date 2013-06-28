#ifndef _KERNELS_H_
#define _KERNELS_H_

#define MAX_KERNEL_DIMENSION %(max_kernel_dimension)d

#include "kernel_base.hpp"

namespace latticeQBP {

%(kernel_source)s

}

#endif
