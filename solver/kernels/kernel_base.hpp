#ifndef _KERNEL_BASE_H_
#define _KERNEL_BASE_H_

////////////////////////////////////////////////////////////////////////////////
// A Kernel wrapper that should be used based on the previous.

#include "array_wrapper.hpp"
#include <algorithm>

template <int d>
static inline Array<int, d> __to_array(std::initializer_list<int> l) {
  Array<int, d> r;

  std::copy(l.begin(), l.end(), r.begin());

  return r;
}


template <int outer_d, int inner_d>
static inline Array<Array<int, inner_d>, outer_d> 
__to_nested_array(std::initializer_list<std::initializer_list<int> > ll) {
  Array<Array<int, inner_d>, outer_d> rr;

  for(auto it = ll.begin(); it != ll.end(); ++it) {
    int i = it - ll.begin();
    std::copy(it->begin(), it->end(), rr[i].begin());
  }

  return rr;
}

template <int _n_dimensions, int _size> class KernelBase {
public:
  static constexpr int kernel_verification = 1;

  static_assert(_size % 2 == 0, "Kernel must be symmetric!");

  typedef Array<int, _n_dimensions> delta_type;
  typedef Array<delta_type, _size> deltas_type;

  inline KernelBase(std::initializer_list<std::initializer_list<int> > _deltas)
    : deltas(__to_nested_array<_size, _n_dimensions>(_deltas))
  {}

  const deltas_type deltas;

  static constexpr size_t n_dimensions = _n_dimensions;
  static constexpr size_t size = _size;
  static constexpr size_t positive_size = _size / 2;

  static inline constexpr size_t kernelSize() { return _size; }
  static inline constexpr size_t kernelPositiveSize() { return _size / 2; }

  static inline constexpr bool isPositiveDirection(unsigned int ei) {
    return ei < kernelPositiveSize(); 
  }

  static inline constexpr int reverseIndex(unsigned int ei) { 
    return ( ( !(_size & (_size - 1)) ) // a multiple of 2
	     ? ((ei + positive_size) & (_size - 1))
	     : (ei + (isPositiveDirection(ei) ? 1 : -1)*kernelPositiveSize()));
    //return isPositiveDirection(ei) ? ei + kernelPositiveSize() : ei - kernelPositiveSize();
  }
};




#endif /* _KERNEL_BASE_H_ */
