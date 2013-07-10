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
    typedef TVRegPathSegment<dtype, TV_PR_Class> RegPathSegment;
    typedef typename RegPathSegment::Mode RPSMode;

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

    vector<RegPathSegment*> node_map_at_lambda_max;

    ////////////////////////////////////////////////////////////////////////////////
    // For managing new regpath instances 

    //. Use deque, as addresses are guaranteed not to change; still
    //have O(1) lookup by index (key)
    deque<RegPathSegment> _regpathsegment_hold;

    RegPathSegment* getNewRegPathSegment() {
      _regpathsegment_hold.emplace_back
        (uint(_regpathsegment_hold.size()), lattice, solver);

      return &_regpathsegment_hold.back();
    }

    RegPathSegment* lookupRPSFromKey(uint key) const {
      return _regpathsegment_hold[key];
    }

    ////////////////////////////////////////////////////////////////////////////////
    // INITIAL: methods for constructing the initial reg paths

    void _constructInitialRegPathsFromSolvedLattice(dtype solved_lamba) {

      if(DEBUG_MODE) {
        for(node_ptr n: lattice) 
          assert(n->keyIsClear());
      }

      // Reserve the 
      node_map_at_lambda_max.resize(lattice.end() - lattice.begin());
      list<RegPathSegment*> active_rps_nodes;

      // Use the lattice for the 
      for(node_ptr seed_node = lattice.begin(); seed_node != lattice.begin(); ++seed_node) {

        RegPathSegment *rps;

        if(unlikely(seed_node->keyIsClear())) {

          dtype lvl = seed_node->level();

          vector<node_ptr> region = solver.walkConnectedRegion
            (seed_node, 
             [lvl](node_ptr nn) { return abs(lvl - nn->level()) <= 1; } );

          rps = getNewRegPathSegment();
        
          rps->setupAsInitial(solved_lamba, region.begin(), region.end());

          assert(rps == lookupRPSFromKey(seed_node->key()));

          active_rps_nodes.push_back(rps);

        } else {
          rps = lookupRPSFromKey(seed_node->key());
          
          // Create the new region; setting all the sections to the
          // correct key. 
        }

        node_map_at_lambda_max[seed_node - lattice.begin()] = rps;
      }

      ////////////////////////////////////////////////////////////////////////////////
      // Now go through and construct all the neighborhood maps.  
      
      for(node_ptr n : lattice) {

        for(uint ei = 0; ei < lattice.kernel_positive_size; ++ei) {
          // blreuhdedon
        }
      }
      
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Building things 
    


    void buildPath

    list<RegPathSegment*> active_rps_nodes;

    

    

  };
}; 

#endif /* _TV_REGPATH_H_ */
