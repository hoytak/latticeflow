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

#include "../common.hpp"
#include "indexing.hpp"

#include <vector>
#include <iostream>

namespace latticeQBP {

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

    LatticeArray(const index_vect& _dimensions, const index_vect& edge_buffers = index_vect(0))
      : indexer( (_dimensions + edge_buffers) + Element<size_t>(0, 3*edge_buffers[0]))
      , bounds(_dimensions)
      , data(indexer.size())
      , begin_ptr(data.begin() + (2*edge_buffers[0])*indexer.stride(0))
      , end_ptr(data.end() - (2*edge_buffers[0])*indexer.stride(0))
    {
    }

    LatticeArray ( LatticeArray && ) = default;

    inline value_ptr begin() const {
      return begin_ptr;
    }

    inline value_ptr end() const {
      return end_ptr;
    }
  
    inline value_ptr ptr(const index_vect& idxv) const {
      return begin() + indexer.getIndex(idxv);
    }

    inline value_ptr ptr(size_t idx) const {
      return begin() + idx;
    }

    template <typename... Pos>
    inline value_ptr operator()(const Pos&... idxargs,
                                typename std::enable_if<sizeof...(Pos) == nd>::type* = 0) {
      return (*this)(index_vect(idxargs...));
    }

    inline value_ptr operator()(const index_vect& idxv) const {
      return ptr(idxv);
    }

    inline value_ptr operator()(const IndexIterator<nd>& it) const {
      return ptr(it.coords());
    }

    inline value_ptr operator()(const BoundedIndexIterator<nd>& it) const {
      return ptr(it.coords());
    }

    inline value_ptr operator()(size_t idx) const {
      return ptr(idx);
    }

    template <typename... Pos>
    inline value_cref operator[](const Pos&... pos) const {
      return *(*this)(pos...);
    }

    template <typename... Pos>
    inline value_ref operator[](const Pos&... pos) {
      return *(*this)(pos...);
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

    template <typename TT> 
    bool withinBounds(const TT* t,
                      typename enable_if<is_base_of<TT, value_type>::value>::type* = 0) const {
      value_direct_cptr vptr = static_cast<value_direct_cptr>(t);
      
      return withinBounds(indexer.getCoords(index(vptr)));
    }

    bool withinBounds(const value_cptr& vptr) const {
      return withinBounds(indexer.getCoords(index(vptr)));
    }

    index_vect shape() const {
      return bounds;
    }

    index_vect fullShape() const {
      return indexer.dimensions();
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

    template <typename TT> 
    bool isValidNode(const TT* t,
                typename enable_if<is_base_of<TT, value_type>::value>::type* = 0) const {
      value_direct_cptr nn = static_cast<value_direct_cptr>(t);
      return (nn >= &(*(data.begin()))) && (nn < &(*(data.end())));
    }

    inline value_ptr resolve(value_direct_ptr vpt) const {
      return begin_ptr + (vpt - &(*begin_ptr));
    }

    inline value_ptr resolve(value_type& vpt) const {
      return resolve(&vpt);
    }

    inline value_ptr resolve(const value_ptr& vpt) const {
      return vpt;
    }

    size_t sizeWithinBounds() const {
      return prod(bounds);
    }

    size_t index(index_vect x) const {
      return indexer.getIndex(x);
    }

    size_t boundedIndex(index_vect x) const {
      return Indexer<nd>(bounds).getIndex(x);
    }

  protected:
    const Indexer<nd> indexer;
    index_vect bounds;
    std::vector<T> data;
    const value_ptr begin_ptr, end_ptr;
  };

};

#ifdef EMACS_FLYMAKE
#include "../kernels/kernels.hpp"
  
  namespace latticeQBP {
    template class LatticeArray<int, 2>;
  };
#endif

#endif /* _LATTICE_H_ */
