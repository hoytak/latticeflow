#ifndef _MN_REDUCTION_POLICIES_H_
#define _MN_REDUCTION_POLICIES_H_

#include "mn_netflow_solver.hpp"
#include "network_flow_node.hpp"
#include "lattice.hpp"
#include "simple_setup.hpp"

template <int _n_dimensions, typename _Kernel, typename _dtype> 
struct _NetFlowReductionMinimizationPolicy {

  // constexpr int n_dimensions() {return _n_dimensions; }

  typedef _Kernel Kernel;
  typedef NetworkFlowNode<_Kernel, _dtype, 2> Node;
  typedef _dtype dtype;

  typedef KernelLattice<Node, _n_dimensions, _Kernel> Lattice;

  typedef typename Node::template NodeFiller<Lattice> Filler;

  typedef NoSetup<Lattice> Setup;
  
  typedef NetflowReductionSolver<dtype, Lattice> Solver;
};


#endif /* _NETWORK_FLOW_POLICY_H_ */
