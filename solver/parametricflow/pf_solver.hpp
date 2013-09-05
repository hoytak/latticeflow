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

    vector<vector<node_ptr> > run() {

      if(HAVE_OUTPUT)
        cout << "Setting up model." << endl;

#if HAVE_OUTPUT
      TimeTracker tt;
      tt.start();
#endif
      // Set up the solver
      NFSolver nf_solver(lattice);
      nf_solver.disablePrinting();

      vector<node_ptr> ls_points;
      typedef typename vector<node_ptr>::iterator vn_iter;

      vector<size_t> boundaries;
      
      ls_points.reserve(lattice.size());

      for(node_ptr n = lattice.begin(); n != lattice.end(); ++n) {
        if(lattice.withinBounds(n)) 
          ls_points.push_back(n);
      }

      struct VNRange {vn_iter start, end;};

      deque<VNRange> queue;

      unsigned int key = 1;

      queue.push_back({ls_points.begin(), ls_points.end()});
      boundaries.push_back(0);

      // keep partitioning the queue
      while(!queue.empty()) {
        auto ls_range = queue.front();
        queue.pop_front();

        assert(ls_range.start != ls_range.end);

        // cout << "running reductions on nodes: ";
        
        // sort(ls_range.start, ls_range.end);

        // for(vn_iter it = ls_range.start; it != ls_range.end; ++it) {
        //   cout << ((*it) - lattice.begin()) << ",";
        // }

        // cout << endl;

        // cerr << "<<<<<<<<<<<<<<<" << '\n';

        // Set them to their initial value
        if(Node::setToMeanReductionPoint(lattice, ls_range.start, ls_range.end)) {
          // cout << "DONE" << endl;
          // cout << "Done with reductions on nodes: ";
          // sort(ls_range.start, ls_range.end);
          // for(vn_iter it = ls_range.start; it != ls_range.end; ++it) {
          //   cout << ((*it) - lattice.begin()) << ",";
          // }
          // cout << endl;
          continue;
        }

        // Run it 
        nf_solver.prepareSection(ls_range.start, ls_range.end, key, false);
        nf_solver.runSection(ls_range.start, ls_range.end, key);

        // Go through and separate out the different levels, if there are any 
        vn_iter it_fwd = ls_range.start;
        vn_iter it_bck = ls_range.end - 1;

        long mid_point = 0;

        assert(it_fwd != it_bck);

        while(true) {
          // cerr << "Fwd: starting; at " << (it_fwd - ls_range.start) << endl;
          // cerr << "Bck starting; at " << (it_bck - ls_range.start) << endl;

          assert(it_fwd != it_bck);

          while( (*it_fwd)->state() == 0 ) {
            ++it_fwd;
            if(unlikely(it_bck == it_fwd))
              goto separation_done;
            // cerr << "ran fwd; at ." << (it_fwd - ls_range.start) << endl;
          }
          
          // cerr << "Done" << endl;

          while( (*it_bck)->state() == 1 ) {
            --it_bck;
            if(unlikely(it_bck == it_fwd)) 
              goto separation_done;
            // cerr << "ran back; at ." << (it_bck - ls_range.start) << endl;
          }
          
          // cerr << "Done" << endl;

          swap(*it_fwd, *it_bck);
          ++it_fwd;

          if(it_fwd == it_bck) {
          separation_done:;
            
            if(it_bck != ls_range.end && (*it_bck)->state() == 0)
              ++it_bck;
            
            if(DEBUG_MODE) {            
              for(auto it = ls_range.start; it != it_bck; ++it) 
                assert((*it)->state() == 0);

              for(auto it = it_bck; it != ls_range.end; ++it) 
                assert((*it)->state() == 1);
            }
  
            mid_point = (it_bck - ls_range.start);
            break;
          }
        }

        // cerr << "Adding mid_point = [" << mid_point 
        //      << "/" << (ls_range.end - ls_range.start) << endl;

        long top_idx = (ls_range.end - ls_range.start); 

        if(mid_point > 0)
          queue.push_back({ls_range.start, ls_range.start + mid_point});

        if(mid_point < top_idx) {
          queue.push_back({ls_range.start + mid_point, ls_range.end});
          
          for(vn_iter it = ls_range.start + mid_point; it != ls_range.end; ++it) 
            (*it)->template flipNode<1>(lattice);
        }

        if(mid_point > 0 && mid_point < top_idx)
          boundaries.push_back(mid_point + (ls_range.start - ls_points.begin()));
        
        nf_solver.cleanupSection(ls_range.start, ls_range.end, key);
        
        // Next key
        ++key;
      }

      // cout << "boundaries = " << boundaries << endl;

      // Go through and reset all the offset stuff
      for(vn_iter it = ls_points.begin(); it != ls_points.end(); ++it) 
        (*it)->setOffset(lattice, 0);

      size_t n_level_sets = boundaries.size();
      boundaries.push_back(ls_points.size());
      sort(boundaries.begin(), boundaries.end());

#if HAVE_OUTPUT
        cout << "Finished running model in " << tt.asString() << "." << endl;
#endif

      // Now just convert these to keysets 
      vector<vector<node_ptr> > level_sets(n_level_sets);
      
      for(size_t i = 0; i < n_level_sets; ++i) {
        level_sets[i] = vector<node_ptr>(ls_points.begin() + boundaries[i],
                                         ls_points.begin() + boundaries[i+1]);
        sort(level_sets[i].begin(), level_sets[i].end());
        
#ifndef NDEBUG
        assert(!level_sets[i].empty());
        for(node_ptr n : level_sets[i])
          assert(lattice.withinBounds(n));
#endif
      }
      
      return level_sets;
    }
  };

  



};

#include "../common/debug_flymake_test.hpp"

#endif /* _MN_NETFLOW_SOLVER_H_ */

