#ifndef _MN_REDUCTION_POLICIES_H_
#define _MN_REDUCTION_POLICIES_H_

#include "pf_solver.hpp"
#include "pf_flow_node.hpp"
#include "../lattices/kernellattice.hpp"
#include "../energy_minimization/energy_base.hpp"
#include "../common.hpp"

namespace latticeQBP {

  template <int _n_dimensions, typename _Kernel, typename _dtype> 
  struct _UnweightedParametricFlowReductionPolicy {

    typedef _Kernel Kernel;
    typedef PFFlowNode<_Kernel, _dtype, PFUnweightedNodePolicy> Node;
    typedef _dtype dtype;

    typedef KernelLattice<Node, _n_dimensions, _Kernel> Lattice;

    typedef typename Node::template NodeFiller<Lattice> Filler;

    typedef NoSetup<Lattice> Setup;
  
    typedef ParametricFlowSolver<dtype, Lattice> Solver;
  };

  // For now, use the same interface as the energy minimization stuff. 

  template <int n_dimensions, typename Kernel, typename dtype = long>
  class LatticeLevelReductions 
    : public LatticeEnergyBase<_UnweightedParametricFlowReductionPolicy<n_dimensions, Kernel, dtype> >
  {

  public:
    typedef LatticeEnergyBase<_UnweightedParametricFlowReductionPolicy<n_dimensions, Kernel, dtype> > Base;
    typedef typename Base::index_vect index_vect;

    LatticeLevelReductions(const index_vect& dimensions)
      : Base(dimensions)
    {}

    double level(const index_vect& node_index) const {
      return Base::lattice(node_index)->level();
    }
  };
};

#endif /* _NETWORK_FLOW_POLICY_H_ */
