#ifndef _TV_PUSH_RELABEL_H_
#define _TV_PUSH_RELABEL_H_

#include "../common.hpp"
#include "../network_flow/push_relabel.hpp"

#include <unordered_set>

namespace latticeQBP {

  using namespace std;  
  
  template <typename dtype, typename KernelLattice>
  class TV_PRFlow : public PRFlow<dtype, KernelLattice, 0> {
  public:
    
    typedef PRFlow<dtype, KernelLattice, 0> Base;
    typedef typename Base::node_ptr node_ptr;
    typedef typename Base::node_cptr node_cptr;
    typedef typename CompType<dtype>::comp_type comp_type;

    typedef unsigned int uint;
    static constexpr uint kernel_size = KernelLattice::kernel_size;    

    TV_PRFlow(KernelLattice& lattice) 
      : Base(lattice)
    {}

    ////////////////////////////////////////////////////////////////////////////////
    // This 

    struct PartitionInfo {
      bool is_on;
      dtype cut_value;
      vector<node_ptr> nodes;
    };

    typedef shared_ptr<PartitionInfo> partitioninfo_ptr;

  public:

    struct CutInfo {
      bool any_cut;
      vector<partitioninfo_ptr> partitions;
      vector<pair<node_ptr, uint> > cut_edges;
    };

    typedef shared_ptr<CutInfo> cutinfo_ptr;

    template <typename NodePtrIterator>  
    cutinfo_ptr runPartitionedSection(const NodePtrIterator& start, 
                                      const NodePtrIterator& end, uint key) {

      if(DEBUG_MODE) {
        for(NodePtrIterator it = start; it != end; ++it) {
          node_ptr n = *it;
          assert(n->matchesKey(key));
          assert(n->height == 0);
          assert(n->state() == 0);
        }
      }

      runSection(start, end, key);

      bool any_on = false;
      for(NodePtrIterator it = start; it != end; ++it) {
        if( (*it)->state() ) {
          any_on = true;
          break;
        }
      }

      cutinfo_ptr cut = cutinfo_ptr(new CutInfo);
      cut->any_cut = any_on;
      
      if(!any_on)
        return cut;

      // First go through and set the state to the proper node.  all
      // these are currently eligible
      for(NodePtrIterator it = start; it != end; ++it) {
        node_ptr n = *it;    

        assert(n->matchesKey(key));
        
        if(!(n->height)) 
          continue;

        // Okay, it's a new partition!!!  Pirate poodles!  

        partitioninfo_ptr pi = partitioninfo_ptr(new PartitionInfo);

        bool is_on = pi->is_on = n->on();
        vector<node_ptr>& nodes = pi->nodes;

        nodes.push_back(n);
        n->height = 1;

        // Go through and fill this partition
        for(size_t n_idx = 0; n_idx < nodes.size(); ++n_idx) {
          node_ptr n2 = nodes[n_idx];
          n2->height = 1;

          for(uint ei = 0; ei < kernel_size; ++ei) {
            node_ptr nn = n2 + Base::step_array[ei];

            if(nn->height != 0 || !nn->matchesKey(key))
              continue;
            
            if(nn->state() == is_on) {
              nodes.push_back(nn);
            } else {

              if(DEBUG_MODE) {
                if(is_on) {
                  assert_lt(n2->level(), nn->level());
                } else {
                  assert_gt(n2->level(), nn->level());
                }
              }

              pi->cut_value += capacityOfSaturated(n2, nn, ei); 

              if(is_on) 
                cut->cut_edges.emplace_back(n2, ei);
            }
          }
        }
        
        cut->partitions.push_back(pi);
      }

      // Finally, clean up 
      for(NodePtrIterator it = start; it != end; ++it) {
        node_ptr n = *it;    

        assert(n->matchesKey(key));
        assert_eq(n->height, 1);

        n->height = 0;
      }
      
      return cut;
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Method for walking the lattice and picking out sections.  This
    // can only be called from an external function; cannot be done while solving

    template <typename UnaryNodeCriterionFunction> 
    vector<node_ptr>  
    walkConnectedRegion(node_ptr seed_node, UnaryNodeCriterionFunction& test_f) {
      
      if(DEBUG_MODE) {
        for(node_ptr n : Base::lattice) {
          assert_eq(n->height, 0);
        }
      }
      
      vector<node_ptr> nodes(1, seed_node);

      size_t visit_index = 0;

      do {
        node_ptr n = nodes[visit_index];

        for(uint ei = 0; ei < kernel_size; ++ei) {

          node_ptr nn = n + Base::step_array[ei];
          
          if(nn->height == 0
             && (pushCapacity(n, nn, ei) + pushCapacity(nn, n, Base::reverseIndex(ei)) > 0)
             && (test_f(n, nn))) {
           
            nn->height = 1;
            nodes.push_back(nn);
          }
        }
      } while( ( (++visit_index) != nodes.size()));

      for(node_ptr n : nodes) 
        n->height = 0;

      if(DEBUG_MODE) {
        for(node_ptr n : Base::lattice) {
          assert_eq(n->height, 0);
        }
      }

      return nodes;
    }

    template <typename ForwardIterator, typename ExtractFunction>
    void getSum(const ForwardIterator& start, const ForwardIterator& end, ExtractFunction& f) {
      comp_type s = 0; 

      for(ForwardIterator it = start; it != end; ++it)
        s += f(*it);

      return s;
    }

    inline bool nodeIsOrphan(node_ptr n) const {
      for(uint ei = 0; ei < Base::lattice.kernel_size; ++ei) {
        node_ptr nn = n + Base::step_array[ei];
          
        if(pushCapacity(n, nn, ei) 
           + pushCapacity(nn, n, Base::reverseIndex(ei)) > 0) {
          return false;
        }

        return true;
      }
    }

    template <typename RegisterNodePair> 
    void constructNeighborhoodKeyPairs(RegisterNodePair& pair_reg) const {

      unordered_set<uint64_t> seen_keys;

      auto gen_key = [](uint key1, uint key2){
        return (uint64_t(min(key1, key2)) << 32) | (uint64_t(max(key1, key2)));
      };

      for(node_ptr n : Base::lattice) {
        for(uint ei = 0; ei < Base::lattice.kernel_positive_size; ++ei) {
          node_ptr nn = n + Base::step_array[ei];

          if(pushCapacity(n, nn, ei) + pushCapacity(nn, n, Base::reverseIndex(ei)) > 0) {
            auto ret = seen_keys.insert(gen_key(n->key(), nn->key()));

            if(ret.second) { 
              // The element was inserted; hasn't been seen before
              pair_reg(n, nn);
            }
          }
        }
      }
    }

    
    template <typename NodePtrIterator>  
    void setRegionToLambda(const NodePtrIterator& start, 
                           const NodePtrIterator& end, 
                           dtype lambda) {
      
      comp_type fv_avg = 0;
      size_t n_sum = 0;

      for(auto it = start; it != end; ++it) {
        node_ptr n = (*it);

        n->setFunctionValue(Base::lattice, 0, lambda); 

        fv_avg += n->fv_predict();
        ++n_sum;
      }

      dtype fv_offset = ceilAverage(fv_avg, n_sum);

      for(auto it = start; it != end; ++it) {
        (*it)->setOffset(Base::lattice, fv_offset);
      }
    }

  };

}; 

#endif /* _TV_PUSH_RELABEL_H_ */
