#ifndef _ENERGY_HPP_
#define _ENERGY_HPP_

#include <iostream>

#include "../common.hpp"
#include "energy_base.hpp"
#include "mn_reduction_policies.hpp"
#include "network_flow_policy.hpp"

namespace latticeQBP {

  using namespace std;

  ////////////////////////////////////////////////////////////////////////////////
  // LatticeEnergyMinimizer

  template <int n_dimensions, typename Kernel, typename dtype = long>
  class LatticeEnergyMinimizer 
    : public LatticeEnergyBase<_SimpleMinimizationPolicy<n_dimensions, Kernel, dtype> >
  {
  public:
    typedef LatticeEnergyBase<_SimpleMinimizationPolicy<n_dimensions, Kernel, dtype> > Base;
    typedef typename Base::index_vect index_vect;

    LatticeEnergyMinimizer(const index_vect& dimensions) : Base(dimensions) {}
  };


  ////////////////////////////////////////////////////////////////////////////////
  // LatticeEnergyAdjuster

  template <int n_dimensions, typename Kernel, typename dtype = long>
  class LatticeEnergyAdjuster
    : public LatticeEnergyBase<_ReductionMinimizationPolicy<n_dimensions, Kernel, dtype> >
  {
  public:
    // Check the types for the simple minimizer
    typedef LatticeEnergyBase<_ReductionMinimizationPolicy<n_dimensions, Kernel, dtype> > Base;
    typedef typename Base::index_vect index_vect;

    LatticeEnergyAdjuster(const index_vect& dimensions) : Base(dimensions) {}
  };

  ////////////////////////////////////////////////////////////////////////////////
  // NetFlowLatticeEnergyReductions

  template <int n_dimensions, typename Kernel, typename dtype = long>
  class LatticeLevelReductions 
    : public LatticeEnergyBase<_NetFlowReductionMinimizationPolicy<n_dimensions, Kernel, dtype> >
  {

  public:
    typedef LatticeEnergyBase<_NetFlowReductionMinimizationPolicy<n_dimensions, Kernel, dtype> > Base;   
    typedef typename Base::index_vect index_vect;

    LatticeLevelReductions(const index_vect& dimensions)
      : Base(dimensions)
    {}

    double level(const index_vect& node_index) const {
      return Base::lattice(node_index)->level();
    }

  };

#ifdef ENABLE_BOYKOV_KOLMOGOROV_GC_CODE 

#include "bk_code_wrapper.hpp"

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





