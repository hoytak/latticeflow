#ifndef _TV_REGPATH_H_
#define _TV_REGPATH_H_

#include "../common.hpp"
#include "tv_flow_node.hpp"
#include "tv_push_relabel.hpp"
#include "tv_regpath_node.hpp"

#include <set>
#include <list>
#include <vector>
#include <deque>
#include <unordered_map>

namespace latticeQBP {
  
  template <typename dtype, typename TV_PR_Class> class TVRegPathManager {
  public:

    typedef typename TV_PR_Class::Lattice Lattice;
    typedef typename TV_PR_Class::node_ptr node_ptr;    
    typedef TVRegPathSegment<dtype, TV_PR_Class> TVRegPathSegment;
    typedef typename TVRegPathSegment::Mode RPSMode;

    TVRegPathManager(Lattice& _lattice, 
                     dtype _max_lambda
                     ) 
      : lattice(_lattice)
      , max_lambda(_max_lambda)
    {
    }

  public:
    void build() {
      
      // We assume that the lattice has been solved with each level
      // being a different partition. 

    }
    
  private:
    ////////////////////////////////////////////////////////////////////////////////
    // Common variables

    Lattice& lattice;
    TV_PR_Class& solver;

    dtype max_lambda;

    vector<TVRegPathSegment*> node_map_at_lambda_max;

    ////////////////////////////////////////////////////////////////////////////////
    // For managing new regpath instances 

    //. Use deque, as addresses are guaranteed not to change; still
    //have O(1) lookup by index (key)
    deque<TVRegPathSegment> _regpathsegment_hold;

    TVRegPathSegment* getNewTVRegPathSegment() {
      _regpathsegment_hold.emplace_back
        (uint(_regpathsegment_hold.size()), lattice, solver);

      return &_regpathsegment_hold.back();
    }

    TVRegPathSegment* lookupRPSFromKey(uint key) const {
      return _regpathsegment_hold[key];
    }

    ////////////////////////////////////////////////////////////////////////////////
    // INITIAL: methods for constructing the initial reg paths

    void _constructInitialRegPathsFromSolvedLattice(dtype solved_lamba) {

      if(DEBUG_MODE) {
        for(node_ptr n: lattice) 
          assert(n->keyIsClear());
      }

      ////////////////////////////////////////////////////////////////////////////////
      // Build an initial set of regions at the primary lambda value

      node_map_at_lambda_max.resize(lattice.end() - lattice.begin());
      vector<TVRegPathSegment*> initial_rps_nodes;

      // The null one is for the dead nodes on the perimeters 
      TVRegPathSegment *null_rps = getNewTVRegPathSegment();

      for(node_ptr seed_node = lattice.begin(); seed_node != lattice.begin(); ++seed_node) {

        TVRegPathSegment *rps;

        if(unlikely(seed_node->keyIsClear())) {

          if(solver.nodeIsOrphan(seed_node)) {
            rps = null_rps;
          } else {

            dtype lvl = seed_node->level();

            vector<node_ptr> region = solver.walkConnectedRegion
              (seed_node, 
               [lvl](node_ptr nn) { return abs(lvl - nn->level()) <= 1; } );

            rps = getNewTVRegPathSegment();
        
            rps->setupAsInitial(solved_lamba, region.begin(), region.end());

            assert(rps == lookupRPSFromKey(seed_node->key()));

            initial_rps_nodes.push_back(rps);
          }

        } else {
          rps = lookupRPSFromKey(seed_node->key());
          
          // Create the new region; setting all the sections to the
          // correct key. 
        }

        node_map_at_lambda_max[seed_node - lattice.begin()] = rps;
      }

      ////////////////////////////////////////////////////////////////////////////////
      // Now go through and construct all the neighborhood maps.  
      auto register_pair = [&](node_ptr n1, node_ptr n2) {
        TVRegPathSegment *rps1 = node_map_at_lambda_max[n1 - lattice.begin()];
        TVRegPathSegment *rps2 = node_map_at_lambda_max[n2 - lattice.begin()];
        
        rps1->neighbors.insert(rps2);
        rps2->neighbors.insert(rps1);
      };
      
      solver.constructNeighborhoodKeyPairs(register_pair);
      

      ////////////////////////////////////////////////////////////////////////////////
      // Now build the entire path

      struct FunPoint {
        dtype lambda;
        enum {Join, Split, SplitUB} mode;
        TVRegPathSegment *rps1, *rps2;
        
        bool operator<(const FunPoint& jp) const {return lambda < jp.lambda;}
      }; 

      priority_queue<FunPoint> run_heap;

      dtype current_lambda = max_lambda;

      ////////////////////////////////////////////////////////////////////////////////
      // Convenience functions to enable quick registration of segments

      auto registerPossibleJoin = 
        [&, run_heap, current_lambda](TVRegPathSegment *rps1, TVRegPathSegment *rps2) {

        dtype join_lambda = TVRegPathSegment::calculateJoins(rps1, rps2, current_lambda);
        if(join_lambda > 0)
          run_heap.push(FunPoint({join_lambda, FunPoint::Join, rps1, rps2}));
      };

      auto registerPossibleSplit = 
        [&, run_heap, current_lambda](TVRegPathSegment *rps) {

        auto sp_info = rps->calculateSplit(current_lambda);
        
        if(sp_info.split_occurs)
          run_heap.push(FunPoint({sp_info.split_lambda, FunPoint::Split, rps, nullptr}));
        else if(sp_info.split_ub != 0)
          run_heap.push(FunPoint({sp_info.split_ub, FunPoint::SplitUB, rps, nullptr}));
      };

      ////////////////////////////////////////////////////////////////////////////////
      // INIT the paths for lookup

      // Init the priority queue for the joins.
      for(TVRegPathSegment* rps1 : initial_rps_nodes) {
        for(TVRegPathSegment* rps2 : rps1->neighbors) {
          if(rps1 < rps2) 
            registerPossibleJoin(rps1, rps2);
        }
      }

      // Now go through and calculate all the splits, putting the
      // split upper bounds and the splits on the heap
      for(TVRegPathSegment* rps : initial_rps_nodes) 
        registerPossibleSplit(rps);

      // They are all in the priority heap at this point
      initial_rps_nodes.clear();

      ////////////////////////////////////////////////////////////////////////////////
      // Now build the reg path
                                                        
      auto activateNewLine = 
        [&, run_heap, current_lambda](TVRegPathSegment* rps) {
        
        assert_equal(current_lambda, rps->rhs_lambda);

        // Check out the joins
        for(TVRegPathSegment* rps2 : rps->neighbors()) 
          registerPossibleJoin(rps, rps2);

        // Calculate potential splits
        registerPossibleSplit(rps);
      };

      auto isStillValid = [&, current_lambda](TVRegPathSegment *rps) {
        assert_neq(rps, nullptr);
        if(rps->lhs_lambda != -1) {
          assert_geq(current_lambda, rps1->lhs_lambda);
          assert(rps1->lhs_mode != TVRegPathSegment::Unset);
          return false;
        }
      };

      // Now we have something for it
      while(!run_heap.empty()) {
        FunPoint fp = run_heap.top();
        run_heap.pop();

        assert_leq(fp.lambda, current_lambda);
        current_lambda = fp.lambda;

        if(!isStillValid(fp.rps1) || (fp.rps2 != nullptr && !isStillValid(fp.rps2)))
          continue;

        switch(fp.mode) {
        case FunPoint::Join:
          TVRegPathSegment* rps1 = fp.rps1;
          TVRegPathSegment* rps2 = fp.rps2;

          TVRegPathSegment* new_rps = getNewTVRegPathSegment();

          TVRegPathSegment::join(current_lambda, new_rps, rps1, rps2);

          rps1->deactivate();
          rps2->deactivate();

          activateNewLine(new_rps);

          break;

        case FunPoint::Split:
          TVRegPathSegment* rps = fp.rps1;
          assert(fp.rps2 == nullptr);

          assert_equal(current_lambda, rps->constructionInfo()->lambda_of_split);

          TVRegPathSegment* new_rps1 = getNewTVRegPathSegment();
          TVRegPathSegment* new_rps2 = getNewTVRegPathSegment();

          rps->applySplit(new_rps1, new_rps2);

          rps->deactivate();

          activateNewLine(new_rps1);
          activateNewLine(new_rps2);

          break;

        case FunPoint::SplitUB:
          TVRegPathSegment* rps = fp.rps1;          
          assert(fp.rps2 == nullptr);

          assert_equal(current_lambda, 
                       rps->constructionInfo()->split_calculation_done_to_lambda);

          registerPossibleSplit(rps);

          break;
        }
      }

      // At this point, everything should be solved...
    }


  };
}; 

#endif /* _TV_REGPATH_H_ */
