#ifndef _ENERGY_BASE_H_
#define _ENERGY_BASE_H_

#include <type_traits>
#include <iostream>
#include "../common.hpp"
#include "../lattices.hpp"

namespace latticeQBP {

  using namespace std;

  template <typename Policy> class LatticeEnergyBase {
  public:
    static_assert(Policy::Kernel::kernel_verification == 1, 
                  "Provided Kernel does not appear to inherit from KernelBase class.");

    static constexpr int n_dimensions = Policy::Kernel::n_dimensions;

    typedef typename Policy::Kernel Kernel;
    typedef Array<size_t, n_dimensions> index_vect;
    typedef typename Policy::dtype dtype;
    typedef typename Policy::Node Node;
    typedef typename Policy::Filler Filler;
    typedef typename Policy::Lattice Lattice;
    typedef typename Policy::Setup Setup;
    typedef typename Policy::Solver Solver;

  protected:

    Lattice lattice;

    LatticeEnergyBase(const index_vect& dimensions) 
      : lattice(dimensions)
    { }

  public:

    ///////////////////////////////////////////////////////////////////////////////
    // Iterator filling; an additional interface to building the lattice

    class UnaryFillingIterator {
    public:

      static constexpr int n_dimensions = Policy::Kernel::n_dimensions;

      inline const index_vect& latticeCoord() const { 
        return node_idx_iter.coords();
      }

      inline size_t nodeIndex() const {
        return node_idx_iter.boundedIndex();
      }

      inline bool done() const { return node_idx_iter.done(); }

      inline void addUnaryPotential(dtype e0, dtype e1) { 
        filler.addE1(lattice.ptr(node_idx_iter.index()), e0, e1);
      }

      inline void addCapacityFromSource(dtype c) {
        filler.addE1(lattice(node_idx_iter.index()), -c);
      }

      inline void addCapacityToSink(dtype c) {
        filler.addE1(lattice(node_idx_iter.index()), c);
      }

      const UnaryFillingIterator& operator++() {
        assert(!node_idx_iter.done());
        node_idx_iter.increment();
        return *this;
      }


    private:
      friend class LatticeEnergyBase;

      UnaryFillingIterator(Lattice& _lattice) 
        : lattice(_lattice)
        , filler(lattice)
        , node_idx_iter(lattice.indexIterator())
      {
      }

      Lattice& lattice;
      Filler filler;
      BoundedIndexIterator<n_dimensions> node_idx_iter;
    };

    class PairwiseFillingIterator {
    public:

      static constexpr int n_dimensions = Policy::Kernel::n_dimensions;

      inline size_t nodeIndexOf1() const {
        return node_idx_iter.boundedIndex();
      }

      inline size_t nodeIndex() const {
        return nodeIndexOf1();
      }

      inline size_t nodeIndexOf2() const {
        return node_idx_target_iters[node_target_iter_index].boundedIndex();
      }

      inline size_t edgeIndex() const {
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

      inline void addUnaryPotentialTo1(dtype e0, dtype e1) { 
        filler.addE1(lattice.ptr(node_idx_iter.index()), e0, e1);
      }

      inline void addUnaryPotentialTo2(dtype e0, dtype e1) { 
        filler.addE1(lattice.ptr(node_idx_target_iters[node_target_iter_index].index()), 
                     e0, e1);
      }

      inline double geocutEdgeWeight() const {
        return lattice.geocutEdgeWeight(edgeIndex());
      }

      inline void addPairwisePotential(dtype e00, dtype e01, dtype e10, dtype e11) {
        filler.addE2(lattice.ptr(node_idx_iter.index()), 
                     lattice.ptr(node_idx_target_iters[node_target_iter_index].index()), 
                     node_target_iter_index,
                     e00, e01, e10, e11);
      }

      inline void addCapacityOfEdge(dtype c) {
        filler.addE2(lattice.ptr(node_idx_iter.index()), 
                     lattice.ptr(node_idx_target_iters[node_target_iter_index].index()), 
                     node_target_iter_index, 
                     0, 0 ,-c);
      }

      inline void addCapacityOfSourceTo1(dtype c) {
        filler.addE1(lattice(node_idx_iter.index()), -c);
      }

      inline void addCapacityOf1ToSink(dtype c) {
        filler.addE1(lattice(node_idx_iter.index()), c);
      }

      inline void addCapacityOfSourceTo2(dtype c) {
        filler.addE1(lattice.ptr(node_idx_target_iters[node_target_iter_index].index()), -c);
      }

      inline void addCapacityOf2ToSink(dtype c) {
        filler.addE1(lattice.ptr(node_idx_target_iters[node_target_iter_index].index()), c);
      }

      const PairwiseFillingIterator& operator++() {
      
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
      friend class LatticeEnergyBase;

      PairwiseFillingIterator(Lattice& _lattice) 
        : lattice(_lattice)
        , filler(lattice)
        , node_idx_iter(lattice.indexIterator())
        , node_idx_target_iters(lattice.fullIndexIterator())
        , node_target_iter_index(0)
      {
        // Advance the non-bounded lattice iterator to the start of the
        // bounded iterator
        for(size_t i = 0; i < lattice.kernelPositiveSize(); ++i)
          node_idx_target_iters[i].jump(lattice.jump(i));
      }

      Lattice& lattice;
      Filler filler;
      BoundedIndexIterator<n_dimensions> node_idx_iter;
      Array<IndexIterator<n_dimensions>, Lattice::kernel_positive_size> node_idx_target_iters;
      size_t node_target_iter_index;
    };

    Lattice& getLattice() {return lattice;}

    const Lattice& getLattice() const {return lattice;}

    PairwiseFillingIterator getPairwiseFillingIterator() {
      return PairwiseFillingIterator(lattice);
    }

    UnaryFillingIterator getUnaryFillingIterator() {
      return UnaryFillingIterator(lattice);
    }

    // Also allow stuff by specifying a single node 
    void addUnaryPotential(const index_vect& coord, dtype e0, dtype e1) {
      Filler(lattice).addE1(lattice.ptr(coord), e0, e1); 
    }

    double run() {
      Setup(lattice).setup();
      Solver solver(lattice);

      TimeTracker tt;
      tt.start();
      solver.run();
      tt.stop();

      return tt.elapsedSeconds();

    }

    bool on(const index_vect& coord) const {
      return lattice(coord)->on();
    }


  };

}

#endif /* _ENERGY_BASE_H_ */
