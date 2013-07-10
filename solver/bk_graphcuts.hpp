#ifndef _BK_CODE_WRAPPER_H_
#define _BK_CODE_WRAPPER_H_

#ifdef ENABLE_BOYKOV_KOLMOGOROV_GC_CODE

#include "lattice.hpp"
#include "bk_graphcuts/bk_energy.hpp"

namespace latticeQBP {

  template <typename Kernel, typename dtype>
  struct _BKGC_Node {
    dtype e1, e0;
    Array<dtype, Kernel::positive_size> e11v, e10v, e01v, e00v;
    bool is_on;
  
    _BKGC_Node() 
      : e1(0), e0(0), e11v(0), e10v(0), e01v(0), e00v(0), is_on(false)
    {}

    inline bool on() const {
      return is_on;
    }
  };

  ////////////////////////////////////////////////////////////////////////////////

  template <class _Lattice, typename dtype> 
  class BKGCFiller {
  public:
    typedef typename _Lattice::value_ptr node_ptr;

    BKGCFiller(_Lattice& _lattice) 
      : lattice(_lattice)
    {}

    inline void addE1(node_ptr n1, dtype _e0, dtype _e1) const {
      n1->e0 = _e0;
      n1->e1 = _e1;
    }

    inline void addE2(node_ptr n1, node_ptr n2, uint ei,
                      dtype e00, dtype e01, dtype e10, dtype e11) const {

      assert(lattice.withinBounds(n1));
      assert(lattice.neighbor(n1, ei) == n2);
      assert(lattice.withinBounds(lattice.neighbor(n1, ei)));
      
      n1->e00v[ei] = e00;
      n1->e01v[ei] = e01;
      n1->e10v[ei] = e10;
      n1->e11v[ei] = e11;
    }

  private:
    _Lattice& lattice;
  };

  ////////////////////////////////////////////////////////////////////////////////

  template <class _KernelLattice, typename dtype>
  class BKGCLattice : public _KernelLattice {
  public:
    typedef typename _KernelLattice::index_vect index_vect;

    BKGCLattice(const index_vect& dimensions) 
      : _KernelLattice(dimensions)
      , bounds(dimensions)
      , bounds_indexer(dimensions)
      , graphcutter(bounds_indexer.size(), _KernelLattice::n_dimensions * bounds_indexer.size())
      , graph_nodes(bounds_indexer.size())
    {
      for(size_t i = 0; i < bounds_indexer.size(); ++i) {
        graph_nodes[i] = graphcutter.add_variable();
      }
    }

    void setup() {

      for(auto it = IndexIterator<_KernelLattice::n_dimensions>(bounds); !it.done(); ++it){
        const auto& n1_coords = it.coords();
        auto n1_ptr = (*this)(n1_coords);
        auto n1_gc_ptr = getGCVar(n1_coords);

        graphcutter.add_term1(n1_gc_ptr, n1_ptr->e0, n1_ptr->e1);

        // cout << "Node: " << n1_coords << ": " << n1_ptr->e0 << ", " << n1_ptr->e1 << endl;

        for(unsigned int i = 0; i < _KernelLattice::kernel_positive_size; ++i) {
          auto n2_coords = this->neighborCoords(n1_ptr, i);

          if(!this->withinBounds(n2_coords))
            continue;

          auto n2_gc_ptr = getGCVar(n2_coords);
	  
          graphcutter.add_term2(n1_gc_ptr, n2_gc_ptr, 
                                n1_ptr->e00v[i], n1_ptr->e01v[i], 
                                n1_ptr->e10v[i], n1_ptr->e11v[i]);
	  
          // cout << "Edge: " << n1_coords << " -- " << n2_coords << ": "
          //      << n1_ptr->e00v[i] << ", " << n1_ptr->e01v[i] << ", " 
          //      << n1_ptr->e10v[i] << ", " << n1_ptr->e11v[i] << endl;
        }
      }
    }

    void run() {
      graphcutter.minimize();

      for(auto it = IndexIterator<_KernelLattice::n_dimensions>(bounds); !it.done(); ++it){
        const auto& n1_coords = it.coords();
        auto n1_ptr = (*this)(n1_coords);
        auto n1_gc_ptr = getGCVar(n1_coords);

        n1_ptr->is_on = (graphcutter.get_var(n1_gc_ptr) != 0);
      }
    }

  private:
    const index_vect bounds;
    Indexer<_KernelLattice::n_dimensions> bounds_indexer;

    typedef Energy<dtype, dtype, dtype> GraphCutter;
    GraphCutter graphcutter;
    vector<typename GraphCutter::Var> graph_nodes;

    typename GraphCutter::Var& getGCVar(const index_vect& idx) {
      return graph_nodes[bounds_indexer[idx]];
    }
  };

  ////////////////////////////////////////////////////////////////////////////////

  template <class _Lattice> class _BKGC_Proxy {
  public:
    _BKGC_Proxy(_Lattice& _lattice) : lattice(_lattice) {}
    void run() { lattice.run(); }
    void setup() { lattice.setup(); }
  private:
    _Lattice& lattice;
  };

  ////////////////////////////////////////////////////////////////////////////////

  template <int _n_dimensions, typename _Kernel, typename _dtype> 
  struct _BKGCMinimizationPolicy {

    typedef _Kernel Kernel;
    typedef _BKGC_Node<Kernel, _dtype> Node;
    typedef _dtype dtype;

    typedef BKGCLattice<KernelLattice<Node, _n_dimensions, Kernel>, _dtype> Lattice;
    typedef BKGCFiller<Lattice, _dtype> Filler;

    typedef _BKGC_Proxy<Lattice> Setup;

    typedef _BKGC_Proxy<Lattice> Solver;
  };

}

#endif /* Enable bk code. */

#endif /* _BK_CODE_WRAPPER_H_ */

