#ifndef _LATTICE_H_
#define _LATTICE_H_

/*********************************************************************************
 *
 *   A generic lattice structure for these operations and
 *   interactions.  
 *
 *   In this case, what we have is two full buffer layers associated
 *   with the first dimension index (this means at, effectively, -1
 *   and n), plus a single, shared buffer layer for all other indices
 *   that wraps around.  It is assumed that all nodes in the buffer
 *   layers are inert; i.e. they don't affect operations elsewhere,
 *   even when included as parts of the optimizations.
 *
 ********************************************************************************/

#include "array_wrapper.hpp"
#include "indexing.hpp"
#include "kernels.hpp"

#include <vector>
#include <iostream>


using namespace std;

class DirDim { 
public:
  DirDim()
    : dim(0), _step(0)
  {}

  DirDim(int dimension,long __step)
    : dim(dimension), _step(__step)
  {
    assert_geq(dimension, 0);
    assert(_step != 0);
  }
  
  int direction() const {
    assert(!isNull());
    return _step > 0 ? 1 : -1;
  }
  
  int dimension() const {
    assert(!isNull());
    return dim;
  }

  long step() const {
    assert(!isNull());
    return _step;
  }

  bool isNull() const {
    return _step == 0;
  } 

  DirDim reversed() const {
    assert(!isNull());
    return DirDim(dim, -_step);
  }
  
private:
  unsigned int dim;
  long _step;
};

////////////////////////////////////////////////////////////////////////////////
// A generic lattice type

template <typename T, int nd> 
class LatticeArray {
public:

private:
  static Array<DirDim, 2*nd> _dirDimArray(const Indexer<nd>& indexer) {
    Array<DirDim, 2*nd> dd_array;
    for(int idx = 0; idx < 2*nd; ++idx) {
      unsigned int dim = indexer.indexDimension(idx);
      int dir = indexer.indexDirection(idx);
      long step = indexer.neighborStep(dim, dir);
      
      dd_array[idx] = DirDim(dim, step);
    }
    
    return dd_array;
  }

public:
  
  typedef Array<size_t, nd> index_vect;
  typedef T value_type;
  typedef typename std::vector<T>::iterator value_ptr;
  typedef typename std::vector<T>::const_iterator value_cptr;

  typedef T* value_direct_ptr;
  typedef const T* value_direct_cptr;
  
  typedef T& value_ref;
  typedef const T& value_cref;

  LatticeArray(const index_vect& _dimensions, const index_vect& edge_buffers = index_vect(1))
    : indexer( (_dimensions + edge_buffers) + Element<size_t>(0, 3*edge_buffers[0]))
    , bounds(_dimensions)
    , data(indexer.size())
    , begin_ptr(data.begin() + (2*edge_buffers[0])*indexer.stride(0))
    , end_ptr(data.end() - (2*edge_buffers[0])*indexer.stride(0))
  {
  }

  value_ptr begin() const {
    return begin_ptr;
  }

  value_ptr end() const {
    return end_ptr;
  }
  
  value_ptr ptr(const index_vect& idxv) const {
    return begin() + indexer.getIndex(idxv);
  }

  value_ptr ptr(size_t idx) const {
    return begin() + idx;
  }

  value_ptr operator()(const index_vect& idxv) const {
    return ptr(idxv);
  }

  value_ptr operator()(size_t idx) const {
    return ptr(idx);
  }

  value_cref operator[](size_t idx) const {
    return *(begin() + idx);
  }
    
  value_ref operator[](size_t idx) {
    return *(begin() + idx);
  }

  size_t index(value_cptr v) const {
    
    assert(v >= begin());
    assert(v < data.end());

    size_t idx = (v - begin());
    assert( begin() + idx == v);

    assert((*this)(indexer.getCoords(idx)) == v);

    return idx;
  }

  size_t index(value_direct_cptr v) const {
    value_direct_cptr start = &(*begin());
    assert(v >= start);

    size_t idx = (v - start);

    assert_lt(idx, size());

    assert( &((*this)[idx]) == v);

    return idx;
  }

  long _tag(value_direct_cptr v) const { 
    return (v - &(*begin()));
  }

  long _tag(value_cptr v) const { 
    return (v - begin());
  }

  size_t size() const {
    return indexer.size() - 2*indexer.stride(0);
  }

  bool isValidNode(value_ptr v) const {
    return (v >= data.begin()) && (v < data.end());
  }
  
  const index_vect& dimensions() const {
    return bounds;
  }

  value_ptr neighbor(size_t index, size_t dim, int step) {
    return data.begin() + indexer.neighborIndex(index, dim, step);
  }

  value_ptr neighbor(const value_ptr& vp, size_t dim, int step) {
    return vp + indexer.neighborStep(dim, step);
  }

  value_cptr neighbor(size_t index, size_t dim, int step) const {
    return data.begin() + indexer.neighborIndex(index, dim, step);
  }

  value_cptr neighbor(const value_cptr& vp, size_t dim, int step) const {
    return vp + indexer.neighborStep(dim, step);
  }

  index_vect getCoords(value_cptr v) const {
    return indexer.getCoords(index(v));
  }

  index_vect getCoords(value_direct_cptr v) const {
    return indexer.getCoords(index(v));
  }

  index_vect getCoords(size_t idx) const {
    return indexer.getCoords(idx);
  }

  int neighborIteratorBound() const {
    return 2*nd;
  }

  unsigned int indexDimension(int idx) const {
    return indexer.indexDimension(idx);
  }

  bool withinBounds(const index_vect& idx) const {

    for(size_t dim = 0; dim < nd; ++dim) 
      if (idx[dim] >= bounds[dim])
	return false;

    return true;
  }

  bool withinBounds(const value_direct_cptr& vptr) const {
    return withinBounds(indexer.getCoords(index(vptr)));
  }

  bool withinBounds(const value_cptr& vptr) const {
    return withinBounds(indexer.getCoords(index(vptr)));
  }

  BoundedIndexIterator<nd> inline indexIterator() const {
    index_vect idx_dim(indexer.dimensions());
    idx_dim[0] = bounds[0];
    return BoundedIndexIterator<nd>(idx_dim, bounds);
  }

  IndexIterator<nd> inline fullIndexIterator() const {
    index_vect idx_dim(indexer.dimensions());
    idx_dim[0] = bounds[0];
    return IndexIterator<nd>(idx_dim);
  }

  bool isValidNode(value_cptr nn) const {
    return (nn >= data.begin()) && (nn < data.end());
  }

  bool isValidNode(size_t idx) const {
    return (idx < size());
  }

  bool isValidNode(value_direct_cptr nn) const {
    return (nn >= &(*(data.begin()))) && (nn < &(*(data.end())));
  }

  inline value_ptr resolve(value_direct_ptr vpt) {
    return begin_ptr + (vpt - &(*begin_ptr));
  }

  inline value_ptr resolve(value_type& vpt) {
    return resolve(&vpt);
  }

  inline value_ptr resolve(const value_ptr& vpt) {
    return vpt;
  }

  size_t sizeWithinBounds() const {
    return prod(bounds);
  }

  
protected:
  const Indexer<nd> indexer;
  index_vect bounds;
  std::vector<T> data;
  const value_ptr begin_ptr, end_ptr;
};

////////////////////////////////////////////////////////////////////////////////
// A more advanced lattice of kernel stuff

template <typename T, int _n_dimensions, typename __Kernel> 
class KernelLattice : public LatticeArray<T, _n_dimensions> {
public:
  typedef __Kernel Kernel;
  typedef LatticeArray<T, _n_dimensions> Base;
  typedef typename Base::index_vect index_vect;
  typedef typename Base::value_ptr value_ptr;
  typedef typename Base::value_cptr value_cptr;
  typedef typename Base::value_direct_ptr value_direct_ptr;
  typedef typename Base::value_direct_cptr value_direct_cptr;

  static constexpr unsigned int n_dimensions = _n_dimensions;

  static constexpr unsigned int kernel_size  = Kernel::size;
  static constexpr unsigned int kernel_positive_size = Kernel::positive_size;

  static constexpr unsigned int kernelSize() { return Kernel::size; }
  static constexpr unsigned int kernelPositiveSize() { return Kernel::positive_size; }

  inline constexpr bool isPositiveDirection(unsigned int ei) {
    return Kernel::isPositiveDirection(ei);
  }

  inline constexpr unsigned int reverseIndex(unsigned int ei) { 
    return Kernel::reverseIndex(ei);
  }

private:
  
  static index_vect _getMaximumJumps() {

    index_vect max_jump(0);
    
    for(const auto& delta : Kernel().deltas) {
      for(size_t i = 0; i < n_dimensions; ++i) {
	max_jump[i] = max(max_jump[i], size_t(abs(delta[i])));
      }
    }

    return max_jump;
  }

public: 

  ////////////////////////////////////////////////////////////////////////////////
  // Now the rest of the lattice stuff

  KernelLattice(const index_vect& _dimensions) 
    : Base(_dimensions, _getMaximumJumps())
  {
    Array<long, kernel_size> jumps;
    
    Kernel kernel;

    // Get the jumps 
    for(size_t i = 0; i < kernel_size; ++i)
      jumps[i] = Base::indexer.neighborStep(kernel.deltas[i]);

    // Sort these to get the 
    array<size_t, kernel_size> index_mapping;
    for(size_t i = 0; i < kernel_size; ++i)
      index_mapping[i] = i;

    // Sort backwards first
    sort(index_mapping.begin(), index_mapping.end(), 
	 [&](size_t i1, size_t i2) { return jumps[i1] > jumps[i2]; } );
    
    // Now the first half should be positive and the upper half
    // positive; reverse the lower so as to get things in increasing order
    reverse(index_mapping.begin(), index_mapping.begin() + kernel_size / 2); 

    for(size_t i = 0; i < kernel_size; ++i) {
      index_jumps[i] = jumps[index_mapping[i]];
      diff_arrays[i] = kernel.deltas[index_mapping[i]];
    }

    // Make sure it's symmetric
    for(size_t i = 0; i < kernel_positive_size; ++i) {
      assert_gt(index_jumps[i], 0);
      assert_equal(index_jumps[i], -index_jumps[i + kernel_positive_size]);
    }
  }
  
  Array<int, n_dimensions> diff(size_t i) const {
    return diff_arrays[i];
  }

  long jump(size_t i) const {
    assert_lt(i, kernel_size);
    return index_jumps[i];
  }

  long positiveJump(size_t i) const {
    assert_lt(i, kernel_positive_size);
    long jump = index_jumps[i];
    assert_gt(jump, 0);
    return jump;
  }

  const Array<long, kernel_size> jumpArray() const {
    return index_jumps;
  }

  size_t neighbor(size_t n_idx, unsigned int ei) const {
    return n_idx + index_jumps[ei];
  }

  value_ptr neighbor(value_ptr n, unsigned int ei) {
    return n + index_jumps[ei];
  }

  value_cptr neighbor(value_cptr n, unsigned int ei) const {
    return n + index_jumps[ei];
  }

  value_direct_ptr neighbor(value_direct_ptr n, unsigned int ei) {
    return n + index_jumps[ei];
  }

  value_direct_cptr neighbor(value_direct_cptr n, unsigned int ei) const {
    return n + index_jumps[ei];
  }

  inline index_vect neighborCoords(value_cptr n, unsigned int ei) const {
    return getCoords(neighbor(n, ei));
  }

  inline index_vect neighborCoords(value_direct_cptr n, unsigned int ei) const {
    return getCoords(neighbor(n, ei));
  }

  inline index_vect neighborCoords(size_t n_idx, unsigned int ei) const {
    return getCoords(neighbor(n_idx, ei));
  }

private:
  
  Array<long, kernel_size> index_jumps;
  typename Kernel::deltas_type diff_arrays;

};

#endif /* _LATTICE_H_ */
