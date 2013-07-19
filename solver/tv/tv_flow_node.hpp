
#ifndef _TV_FLOW_NODE_H_
#define _TV_FLOW_NODE_H_

#include "../common.hpp"
#include "../parametricflow/pf_flow_node.hpp"

#include <algorithm>

namespace latticeQBP {

  using namespace std;

  typedef PFScaledUnweightedNodePolicy TVUnweightedNodePolicy;

  // At this point, we 

  template <class Kernel, typename dtype> class TVFlowNode 
    : public PFFlowNode<Kernel, dtype, TVUnweightedNodePolicy>
  {
  };

};

#ifdef EMACS_FLYMAKE

#include "../kernels/kernels.hpp"

namespace latticeQBP {
  template class TVFlowNode<Star2d_4, long>;
};

#endif

#endif /* _TV_FLOW_NODE_H_ */
