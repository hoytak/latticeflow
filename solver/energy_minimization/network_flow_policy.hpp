#ifndef _NETWORK_FLOW_POLICY_H_
#define _NETWORK_FLOW_POLICY_H_

#include "../network_flow/push_relabel.hpp" 
#include "../network_flow/network_flow_node.hpp"
#include "../lattices/kernellattice.hpp"
#include "../common.hpp"
#include "../bk_graphcuts.hpp"

namespace latticeQBP {

  template <int _n_dimensions, typename _Kernel, typename _dtype> 
  struct _SimpleMinimizationPolicy {

    // constexpr int n_dimensions() {return _n_dimensions; }

    typedef _Kernel Kernel;
    typedef NetworkFlowNode<Kernel, _dtype, 0> Node;
    typedef _dtype dtype;

    typedef KernelLattice<Node, _n_dimensions, Kernel> Lattice;
    typedef typename Node::template NodeFiller<Lattice> Filler;

    typedef NoSetup<Lattice> Setup;
  
    typedef PRFlow<dtype, Lattice, 0> Solver;
  };

  template <int _n_dimensions, typename _Kernel, typename _dtype> 
  struct _ReductionMinimizationPolicy {

    // constexpr int n_dimensions() {return _n_dimensions; }

    typedef _Kernel Kernel;
    typedef NetworkFlowNode<Kernel, _dtype, 0> Node;
    typedef _dtype dtype;

    typedef KernelLattice<Node, _n_dimensions, Kernel> Lattice;

    typedef typename Node::template NodeFiller<Lattice> Filler;

    typedef NoSetup<Lattice> Setup;
  
    typedef PRFlow<dtype, Lattice, 0> Solver;
  };
}

#endif /* _NETWORK_FLOW_POLICY_H_ */
