
#ifndef _ENERGY_HPP_
#define _ENERGY_HPP_

#include <iostream>

#include "energy_base.hpp"
#include "bk_graphcuts.hpp"
#include "../network_flow/push_relabel.hpp" 
#include "../network_flow/network_flow_node.hpp"
#include "../lattices/kernellattice.hpp"
#include "../common.hpp"

namespace latticeQBP {

  using namespace std;

  template <int _n_dimensions, typename _Kernel, typename _dtype> 
  struct _ReductionMinimizationPolicy {

    // constexpr int n_dimensions() {return _n_dimensions; }

    typedef _Kernel Kernel;
    typedef NetworkFlowNode<Kernel, _dtype, NFNodeDefaultPolicy> Node;
    typedef _dtype dtype;

    typedef KernelLattice<Node, _n_dimensions, Kernel> Lattice;

    typedef typename Node::template NodeFiller<Lattice> Filler;

    typedef NoSetup<Lattice> Setup;
  
    typedef PRFlow<dtype, Lattice, 0> Solver;
  };

  ////////////////////////////////////////////////////////////////////////////////
  // LatticeEnergyAdjuster

  template <int n_dimensions, typename Kernel, typename dtype = long>
  class LatticeEnergyMinimizer
    : public LatticeEnergyBase<_ReductionMinimizationPolicy<n_dimensions, Kernel, dtype> >
  {
  public:
    // Check the types for the simple minimizer
    typedef LatticeEnergyBase<_ReductionMinimizationPolicy<n_dimensions, Kernel, dtype> > Base;
    typedef typename Base::index_vect index_vect;

    LatticeEnergyMinimizer(const index_vect& dimensions) : Base(dimensions) {}
  };

  ////////////////////////////////////////////////////////////////////////////////
  // NetFlowLatticeEnergyReductions

#ifdef ENABLE_BOYKOV_KOLMOGOROV_GC_CODE 

  template <int n_dimensions, typename Kernel, typename dtype = long>
  class BKGCEnergyMinimizer 
    : public LatticeEnergyBase<_BKGCMinimizationPolicy<n_dimensions, Kernel, dtype> >
  {
  public:
    typedef LatticeEnergyBase<_BKGCMinimizationPolicy<n_dimensions, Kernel, dtype> > Base;
    typedef typename Base::index_vect index_vect;

    BKGCEnergyMinimizer(const index_vect& dimensions)
      : Base(dimensions)
    {}
  };

#endif

}

#endif /* _ENERGY_H_ */





