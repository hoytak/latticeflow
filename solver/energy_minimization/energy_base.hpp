#ifndef _ENERGY_BASE_H_
#define _ENERGY_BASE_H_

#include <type_traits>
#include <iostream>
#include "../common.hpp"
#include "../lattices/kernellattice.hpp"

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
    typedef typename Lattice::VertexIterator VertexIterator;
    typedef typename Lattice::EdgeIterator EdgeIterator;

  protected:

    Lattice lattice;

    LatticeEnergyBase(const index_vect& dimensions) 
      : lattice(dimensions)
    { }

  public:

    ///////////////////////////////////////////////////////////////////////////////
    // Iterator filling; an additional interface to building the lattice

    class UnaryFillingIterator : public VertexIterator {
    public:

      inline void addUnaryPotential(dtype e0, dtype e1) { 
        filler.addE1(Base::node(), e0, e1);
      }

      inline void addCapacityFromSource(dtype c) {
        filler.addE1(Base::node(), -c);
      }

      inline void addCapacityToSink(dtype c) {
        filler.addE1(Base::node(), c);
      }

    private:
      friend class LatticeEnergyBase;
      typedef VertexIterator Base;

      UnaryFillingIterator(Lattice& _lattice) 
        : VertexIterator(_lattice)
        , filler(lattice)
      {
      }

      Filler filler;
    };

    class PairwiseFillingIterator : public EdgeIterator {
    public:

      inline void addUnaryPotentialTo1(dtype e0, dtype e1) { 
        filler.addE1(Base::node1(), e0, e1);
      }

      inline void addUnaryPotentialTo2(dtype e0, dtype e1) { 
        filler.addE1(Base::node2(), e0, e1);
      }

      inline void addPairwisePotential(dtype e00, dtype e01, dtype e10, dtype e11) {
        filler.addE2(Base::node1(), Base::node2(), Base::edgeIndex(), e00, e01, e10, e11);
      }

      inline void addCapacityOfEdge(dtype c) {
        filler.addE2(Base::node1(), Base::node2(), Base::edgeIndex(), 0, 0 ,-c);
      }

      inline void addCapacityOfSourceTo1(dtype c) {
        filler.addE1(Base::node1(), -c);
      }

      inline void addCapacityOf1ToSink(dtype c) {
        filler.addE1(Base::node1(), c);
      }

      inline void addCapacityOfSourceTo2(dtype c) {
        filler.addE1(Base::node2(), -c);
      }

      inline void addCapacityOf2ToSink(dtype c) {
        filler.addE1(Base::node2(), c);
      }

    private:
      friend class LatticeEnergyBase;
      typedef EdgeIterator Base;

      PairwiseFillingIterator(Lattice& _lattice) 
        : Base(_lattice)
        , filler(_lattice)
      {
      }

      Filler filler;
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
