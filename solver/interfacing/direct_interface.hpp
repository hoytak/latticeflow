#ifndef _DIRECT_INTERFACE_H_
#define _DIRECT_INTERFACE_H_

#include "../interface.hpp"
#include "../common.hpp"

////////////////////////////////////////////////////////////////////////////////
// Now the interface high level class.
namespace latticeQBP {

  template <class _Solver>
  class LatticeFlowDirectInterface : public LatticeFlowInterface {
  public:
    typedef _Solver Solver;

    static constexpr size_t n_dimensions = Solver::n_dimensions;

    typedef typename Solver::Kernel Kernel;
    typedef typename Solver::index_vect index_vect;
    typedef typename Solver::Lattice Lattice;
    typedef typename Solver::Filler Filler;

  private:
    Solver solver;
    double run_time; 

  public:

    LatticeFlowDirectInterface(index_vect dimensions)
      : solver(dimensions.begin()), run_time(0)
    {
      const size_t n_filling_edges = Lattice::kernelPositiveSize();
    
      vector<int> edge_deltas(n_dimensions * n_filling_edges);
    
      for(size_t j = 0; j < n_filling_edges; ++j) {
        const auto& diff = solver.getLattice().diff(j);

        for(size_t d = 0; d < n_dimensions; ++d) {
          edge_deltas[j*n_dimensions + d] = diff[d];
        }
      }

      vector<size_t> dim_vect(dimensions.begin(), dimensions.end());

      LatticeFlowInterface::init(n_dimensions, n_filling_edges, dim_vect, edge_deltas);
    }

    void run() {
      TimeTracker tt;
      tt.start();

      solver.run();

      tt.stop();

      run_time = tt.elapsedSeconds();
    }

    double runTime() {
      return run_time;
    }


    void addEnergyPotentials(long *unary_energy_array, long *binary_energy_array) {
    
      // Here we have the filler for the 
      for(auto fill_iter = solver.getUnaryFillingIterator(); 
          !fill_iter.done(); ++fill_iter) {

        size_t node_index = fill_iter.nodeIndex();

        assert_lt(node_index, nNodes());

        fill_iter.addUnaryPotential(unary_energy_array[2*node_index + 0],
                                    unary_energy_array[2*node_index + 1]);
      }

      for(auto fill_iter = solver.getPairwiseFillingIterator();
          !fill_iter.done(); ++fill_iter) {
      
        size_t node_index = fill_iter.nodeIndex();
        size_t edge_index = fill_iter.edgeIndex();
	
        assert_lt(node_index, nNodes());
        assert_lt(edge_index, nFillingEdges());

        size_t base_index = node_index * (4*Lattice::kernelPositiveSize()) + 4*edge_index;

        fill_iter.addPairwisePotential(binary_energy_array[base_index + 0],
                                       binary_energy_array[base_index + 1],
                                       binary_energy_array[base_index + 2],
                                       binary_energy_array[base_index + 3]);
      }
    }

    vector<int> getCut() const {
      const Lattice& lattice = solver.getLattice();

      vector<int> cut(lattice.sizeWithinBounds());

      for(auto node_iter = lattice.indexIterator(); !node_iter.done(); ++node_iter)
        cut[node_iter.boundedIndex()] = lattice[node_iter.index()].on() ? 1 : 0;

      return cut;
    }

  };

}


#endif /* _DIRECT_INTERFACE_H_ */
