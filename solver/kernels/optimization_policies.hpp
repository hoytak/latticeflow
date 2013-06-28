#ifndef _KERNEL_OPTIMIZATION_POLICES_H_
#define _KERNEL_OPTIMIZATION_POLICES_H_

template <class Kernel> struct KernelOptimizationPolicy {

  // Initialization
  static constexpr bool init_hotstart()       { return true; }
  static constexpr bool init_quickflow()      { return true; }
  static constexpr bool init_rerun_hotstart() { return false; }

  static constexpr bool use_tiny_kernel_mode()
  {
    return (Kernel::size > 48);
  }

  static constexpr bool use_fixed_sized_sorter() {
    return (Kernel::size <= 48);
  }

  static constexpr bool use_truncated_sorting() {
    return (Kernel::size > 64);
  }

  // Running alternative modes
  static constexpr bool run_topology_restructure() { return false; }

  static constexpr size_t topology_restructure_trigger_threshold() 
  { 
    return 20; 
  }
  
  static constexpr size_t topology_restructure_sink_fill_bonus() 
  { 
    return 20; 
  }

  // Set to 100 to disable
  static constexpr size_t quickflow_percent_after_topology_restructure_threshold()
  { 
    return 10;
  }

};

#include "tuned_optimization_policies.hpp"

#endif /* _KERNEL_OPTIMIZATION_POLICES_H_ */
