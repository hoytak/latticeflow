#ifndef _KERNEL_LATTICE_H_
#define _KERNEL_LATTICE_H_

#include "lattice.hpp"

////////////////////////////////////////////////////////////////////////////////
// A more advanced lattice of kernel stuff

namespace latticeQBP {

  using namespace std;

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
    typedef typename Base::value_type value_type;

    static constexpr unsigned int n_dimensions = _n_dimensions;

    static constexpr unsigned int kernel_size  = Kernel::size;
    static constexpr unsigned int kernel_positive_size = Kernel::positive_size;

    static constexpr unsigned int kernelSize() { return Kernel::size; }
    static constexpr unsigned int kernelPositiveSize() { return Kernel::positive_size; }

    static constexpr bool isPositiveDirection(unsigned int ei) {
      return Kernel::isPositiveDirection(ei);
    }

    static constexpr unsigned int reverseIndex(unsigned int ei) { 
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
        geocut_edge_weights[i] = kernel.geocut_edge_weights[index_mapping[i]];
      }

      // Make sure it's symmetric
      for(size_t i = 0; i < kernel_positive_size; ++i) {
        assert_gt(index_jumps[i], 0);
        assert_equal(index_jumps[i], -index_jumps[i + kernel_positive_size]);
      }
    }
    
    KernelLattice(KernelLattice&& kl) 
      : Base(kl)
      , index_jumps(kl.index_jumps)
      , geocut_edge_weights(kl.geocut_edge_weights)
      , diff_arrays(kl.diff_arrays)
    {}
  
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

    value_cptr neighbor(value_cptr n, unsigned int ei) const {
      return n + index_jumps[ei];
    }

    value_ptr neighbor(value_ptr n, unsigned int ei) {
      return n + index_jumps[ei];
    }


    template <typename TT> 
    TT* neighbor(TT* vptr, unsigned int ei,
                typename enable_if<is_base_of<TT, value_type>::value 
                && !is_const<TT>::value>::type* = 0)
    {
      return static_cast<value_direct_ptr>(vptr) + index_jumps[ei];
    }

    template <typename TT> 
    const TT* neighbor(const TT* vptr, unsigned int ei,
                      typename enable_if<is_base_of<TT, value_type>::value>::type* = 0) const
    {
      return static_cast<value_direct_cptr>(vptr) + index_jumps[ei];
    }

    inline index_vect neighborCoords(value_direct_cptr n, unsigned int ei) const {
      return this->getCoords(neighbor(n, ei));
    }

    inline index_vect neighborCoords(value_cptr n, unsigned int ei) const {
      return this->getCoords(neighbor(n, ei));
    }

    inline index_vect neighborCoords(size_t n_idx, unsigned int ei) const {
      return this->getCoords(neighbor(n_idx, ei));
    }

    inline double neighborL2Dist(unsigned int ei) const {
      index_vect src = this->getCoords(this->begin());
      index_vect dest = this->neighborCoords(this->begin(), ei);
      long r = 0;

      for(unsigned int i = 0; i < n_dimensions; ++i) {
        long v = long(src[i]) - long(dest[i]);
        r += v*v;
      }

      return sqrt(double(r));
    }

    inline double geocutEdgeWeight(size_t i) const {
      return geocut_edge_weights[i];
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Iteration

    class VertexIterator {
    public:
      typedef value_ptr node_ptr;    
    
      VertexIterator(KernelLattice& _lattice) 
        : lattice(_lattice)
        , node_idx_iter(lattice.indexIterator())
      {}

      inline node_ptr node() const {
        return lattice(node_idx_iter.index());
      }

      inline const index_vect& latticeCoord() const { 
        return node_idx_iter.coords();
      }

      inline size_t nodeIndex() const {
        return node_idx_iter.boundedIndex();
      }

      inline size_t internalIndex() const {
        return node_idx_iter.index();
      }

      inline bool done() const { return node_idx_iter.done(); }

      VertexIterator& operator++() {
        assert(!node_idx_iter.done());
        node_idx_iter.increment();
        return *this;
      }

    private:
      KernelLattice& lattice;
      BoundedIndexIterator<n_dimensions> node_idx_iter;
    };

    class EdgeIterator {
    public:

      typedef  value_ptr node_ptr;    

      EdgeIterator(KernelLattice& _lattice) 
        : lattice(_lattice)
        , node_idx_iter(lattice.indexIterator())
        , node_idx_target_iters(lattice.fullIndexIterator())
        , node_target_iter_index(0)
      {
        // Advance the non-bounded lattice iterator to the start of the
        // bounded iterator
        for(size_t i = 0; i < lattice.kernelPositiveSize(); ++i)
          node_idx_target_iters[i].jump(lattice.jump(i));
      }

      inline size_t nodeIndexOf1() const {
        return node_idx_iter.boundedIndex();
      }

      inline size_t nodeIndex() const {
        return nodeIndexOf1();
      }

      inline size_t nodeIndexOf2() const {
        return lattice.Base::boundedIndex(latticeCoordOf2());
      }

      node_ptr node1() const {
        return lattice.ptr(node_idx_iter.index());
      }

      node_ptr node2() const {
        return lattice.ptr(node_idx_target_iters[node_target_iter_index].index());
      }

      inline unsigned int edgeIndex() const {
        return node_target_iter_index;
      }

      inline const index_vect& latticeCoordOf1() const { 
        return node_idx_iter.coords();
      }

      inline const index_vect& latticeCoordOf2() const { 
        return node_idx_target_iters[node_target_iter_index].coords(); 
      }
    
      inline double L2_SquaredDistance() const { 
        return dist2(latticeCoordOf1(), latticeCoordOf2());
      }

      inline double L2_Distance() const { 
        return sqrt(L2_SquaredDistance());
      }
    
      inline bool done() const { return node_idx_iter.done(); }

    
      inline double geocutEdgeWeight() const {
        return lattice.geocutEdgeWeight(edgeIndex());
      }

      const EdgeIterator& operator++() {
      
        do {
          if(unlikely(done()))
            return *this;

          ++node_target_iter_index;

          if(node_target_iter_index >= lattice.kernel_positive_size) {

            assert(!node_idx_iter.done());

            node_target_iter_index = 0;

            // Increment everything
            size_t step = node_idx_iter.increment();
	  
            if(step == 1) {
              for(auto& it : node_idx_target_iters)
                ++it;
            } else {
              for(auto& it : node_idx_target_iters) {
                it.jump(step);
              }
            }
          }

        } while(!lattice.withinBounds(node_idx_target_iters[node_target_iter_index].coords()));

        return *this;
      }

    private:
      KernelLattice& lattice;
      BoundedIndexIterator<n_dimensions> node_idx_iter;
      Array<IndexIterator<n_dimensions>, KernelLattice::kernel_positive_size> node_idx_target_iters;
      size_t node_target_iter_index;
    };

    VertexIterator vertexIterator() {
      return VertexIterator(*this);
    }

    EdgeIterator edgeIterator() {
      return EdgeIterator(*this);
    }

  private:
  
    Array<long, kernel_size> index_jumps;
    Array<double, kernel_size> geocut_edge_weights;
    typename Kernel::deltas_type diff_arrays;

  };
};

#endif /* _KERNEL_LATTICE_H_ */
