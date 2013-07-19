#ifndef _MN_NETFLOW_SOLVER_H_
#define _MN_NETFLOW_SOLVER_H_

#include "../network_flow/push_relabel.hpp"
#include <list>
#include <deque>
#include "../common.hpp"

namespace latticeQBP {

  using std::list;
  using std::deque;

  template <typename dtype, typename KernelLattice> 
  class ParametricFlowSolver {
  private:
    KernelLattice& lattice;
    typedef PRFlow<dtype, KernelLattice, 0> NFSolver;
    
    typedef typename KernelLattice::value_ptr  node_ptr;
    typedef typename KernelLattice::value_type Node;

    static constexpr uint kernel_size = KernelLattice::kernel_size;
    static constexpr uint kernel_positive_size = KernelLattice::kernel_positive_size;

  public:

    ParametricFlowSolver(KernelLattice& _lattice)
      : lattice(_lattice)
    {
    }

    void run() {

      if(HAVE_OUTPUT)
        cout << "Setting up model." << endl;

#if HAVE_OUTPUT
      TimeTracker tt;
      tt.start();
#endif
      // Set up the solver
      NFSolver nf_solver(lattice);

      vector<node_ptr> ls_points;
      typedef typename vector<node_ptr>::iterator vn_iter;

      vector<size_t> boundaries(lattice.size()); 
      
      ls_points.reserve(lattice.size());

      for(node_ptr n = lattice.begin(); n != lattice.end(); ++n) {
        if(lattice.withinBounds(n)) 
          ls_points.push_back(n);
      }

      struct VNRange {vn_iter start, end;};

      deque<VNRange> queue;

      unsigned int key = 0;

      queue.push_back({ls_points.begin(), ls_points.end()});
      boundaries.push_back(0);

      // keep partitioning the queue
      while(!queue.empty()) {
        auto ls_range = queue.front();
        queue.pop_front();

        assert(ls_range.start != ls_range.end);

        // Set them to their initial value
        if(Node::setToMeanReductionPoint(lattice, ls_range.start, ls_range.end))
          continue;

        // Run it 
        nf_solver.prepareSection(ls_range.start, ls_range.end, key, false);
        nf_solver.runSection(ls_range.start, ls_range.end, key);

        // Go through and separate out the different levels, if there are any 
        vn_iter it_fwd = ls_range.start;
        vn_iter it_bck = ls_range.end - 1;

        long mid_point = 0;

        while(true) {
          while( (*it_fwd)->state() == 0 ) {
            ++it_fwd;
            if(unlikely(it_bck == it_fwd)) {
              goto separation_done;
            }
          }

          while( (*it_bck)->state() == 1 ) {
            --it_bck;
            if(unlikely(it_bck == it_fwd)) 
              goto separation_done;
          }
          
          swap(*it_fwd, *it_bck);
          ++it_fwd;
          if(it_fwd == it_bck) {
          separation_done:;
            
            if((*it_bck)->state() == 0)
              ++it_bck;
            
            if(DEBUG_MODE) {            
              for(auto it = ls_range.start; it != it_bck; ++it) 
                assert((*it)->state() == 0);

              for(auto it = it_bck; it != ls_range.end; ++it) 
                assert((*it)->state() == 1);
            }
  
            mid_point = (it_fwd - ls_range.start);
            break;
          }
        }

        if(mid_point != 0 && mid_point != (ls_range.end - ls_range.start)) {
          queue.push_back({ls_range.start, ls_range.start + mid_point});
          queue.push_back({ls_range.start + mid_point, ls_range.end});
        }

        // Now clean up the graph
        nf_solver.cleanupSection(ls_range.start, ls_range.end, key);
        
        // Next key
        ++key;
      }

      if(HAVE_OUTPUT)
        cout << "Finished running model in " << tt.asString() << "." << endl;
    }
  };
};

#ifdef EMACS_FLYMAKE

#include "../kernels/kernels.hpp"
#include "../lattices/kernellattice.hpp"
#include "pf_flow_node.hpp"

namespace latticeQBP {
  typedef KernelLattice<PFFlowNode<Star2d_4, long, PFUnweightedNodePolicy>, 2, Star2d_4> _TV_PRFlow_TestLattice;
  template class ParametricFlowSolver<long, _TV_PRFlow_TestLattice>;
};

#endif

#endif /* _MN_NETFLOW_SOLVER_H_ */
