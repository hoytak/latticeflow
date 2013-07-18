#ifndef _MN_NETFLOW_SOLVER_H_
#define _MN_NETFLOW_SOLVER_H_

#include "../network_flow/push_relabel.hpp"
#include "../network_flow/network_flow_node.hpp"
#include <list>
#include <deque>
#include "../common.hpp"

namespace latticeQBP {

  using std::list;
  using std::deque;

  template <typename dtype, typename KernelLattice> 
  class NetflowReductionSolver {
  private:
    KernelLattice& lattice;
    typedef typename KernelLattice::value_ptr node_ptr;

    static constexpr uint kernel_size = KernelLattice::kernel_size;
    static constexpr uint kernel_positive_size = KernelLattice::kernel_positive_size;

  public:

    NetflowReductionSolver(KernelLattice& _lattice)
      : lattice(_lattice)
    {
    }

    void run() {

      if(HAVE_OUTPUT)
        cout << "Setting up model." << endl;

      TimeTracker tt;
      tt.start();

      PRFlow<dtype, KernelLattice, 0> pr_flow(lattice);
      vector<dtype> reduction_adjustments(lattice.size(), 0);
      deque<list<node_ptr> > queue;

      // Run this by successive zero-runs in order to keep the number
      // of flipping flips down

      auto do_zero_drill = [&reduction_adjustments, &queue, &pr_flow, this]
        (list<node_ptr>& nl)
        {
          while(true) {
	   
            // Set the correct state
            if(nl.empty()) 
              return;

            dtype r_total = 0;
            dtype r_min = nl.front()->reduction;
            dtype r_max = nl.front()->reduction;

            for(node_ptr n : nl) {
              assert_equal(n->state(), 0);
              r_total += n->reduction;
              r_min = min(r_min, n->reduction);
              r_max = max(r_max, n->reduction);
            }

            // Now, go through and redo it

            dtype r_shift = dtype(ceil(double(r_total) / nl.size()));
	   
            // need to prevent numerical issues here

            if(r_max - r_min <= 1) {
              // We're done

              for(node_ptr n : nl)
                n->template flipNode<0>(lattice);

              return;
            }

            // b/c of ceil, this is enough
            if(r_shift == r_max)
              r_shift = r_max - 1;
	     
            for(node_ptr n : nl) {
              n->adjustReduction(-r_shift);
              reduction_adjustments[lattice.index(n)] -= r_shift;
            }

            pr_flow.prepareSection(nl.begin(), nl.end());

            pr_flow.runSection(nl.begin(), nl.end());

            if(likely(!queue.back().empty()))
              queue.push_back(list<node_ptr>());

            for(auto it = nl.begin(); it != nl.end();) {
              node_ptr n = *it;

              if(n->state()) {
                queue.back().splice(queue.back().end(), nl, it++);
              } else
                ++it;
            }
          }
        };
     
      // Now, do a first run with nothing fancy at all
      {
        pr_flow.run();

        list<node_ptr> zero_list;

        queue.push_back(list<node_ptr>());
       
        for(node_ptr n = lattice.begin(); n != lattice.end(); ++n) {
          if(lattice.withinBounds(n)) {
            if(!(n->state()))
              zero_list.push_back(n);
            else
              queue.back().push_back(n);
          }
        }

        do_zero_drill(zero_list);
      }

      // So, at this point, the queue should be ready to go
      while(!queue.empty()) {
       
        for(node_ptr n : queue.front()) {
          assert(n->state());
          n->template flipNode<1>(lattice);
        }
       
        do_zero_drill(queue.front());
       
        queue.pop_front();
      }

      for(auto& n : lattice) {
        if(lattice.withinBounds(&n)) {
          assert_equal(n.state(), 1);
          n.adjustReduction(-reduction_adjustments[lattice.index(&n)]);
        }
      }

      // Now, just go back thru and make sure that all reductions are accounted for
#ifndef NDEBUG
      for(auto& n : lattice) {
        assert_equal(n.reduction_shift, 0);
      }
#endif

      if(HAVE_OUTPUT)
        cout << "Finished running model in " << tt.asString() << "." << endl;

    }
  };
};
#endif /* _MN_NETFLOW_SOLVER_H_ */
