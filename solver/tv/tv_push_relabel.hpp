#ifndef _TV_PUSH_RELABEL_H_
#define _TV_PUSH_RELABEL_H_

#include "../common.hpp"
#include "../network_flow/push_relabel.hpp"

#include <unordered_set>

namespace latticeQBP {

  using namespace std;  

  template <typename dtype, typename _KernelLattice>
  class TV_PRFlow : public PRFlow<dtype, _KernelLattice, 0> {
  public:
    
    typedef PRFlow<dtype, _KernelLattice, 0> Base;
    typedef typename Base::Lattice Lattice;
    typedef typename Base::node_ptr node_ptr;
    typedef typename Base::node_cptr node_cptr;
    typedef typename Base::Node Node;
    typedef typename CompType<dtype>::Type comp_type;

    typedef unsigned int uint;
    static constexpr uint kernel_size = Lattice::kernel_size;

    TV_PRFlow(Lattice& lattice) 
      : Base(lattice)
    {}

    ////////////////////////////////////////////////////////////////////////////////
    // This 

    struct PartitionInfo {
      bool is_on;
      vector<node_ptr> nodes;
    };

    typedef shared_ptr<PartitionInfo> partitioninfo_ptr;

  public:

    struct CutInfo {
      CutInfo() : partitions({nullptr, nullptr}) {}

      bool any_cut;
      Array<partitioninfo_ptr, 2> partitions;
      dtype cut_value;
      vector<pair<node_ptr, uint> > cut_edges;
      Array<dtype, 2> qii_sum, gamma_sum;
    };

    typedef shared_ptr<CutInfo> cutinfo_ptr;

#if !defined(NDEBUG) && ENABLE_EXPENSIVE_CHECKS

    template <typename NodePtrIterator> 
    void checkPartitionedSection(const NodePtrIterator& start, 
                                 const NodePtrIterator& end, uint key) {


      for(NodePtrIterator it = start; it != end; ++it) {
        node_ptr n = *it;
        assert(n->matchesKey(key));
        assert(n->height == 0);
        assert(n->state() == 0);
      }

      set<node_ptr> node_check_set(start, end);

      for(node_ptr n = Base::lattice.begin(); n != Base::lattice.end(); ++n) {
        if(node_check_set.find(n) == node_check_set.end()) {
          assert(!n->matchesKey(key));
        }
      }
    }
#else 
    template <typename NodePtrIterator> 
    static inline void checkPartitionedSection(const NodePtrIterator&, 
                                               const NodePtrIterator&, uint)
    {}
#endif

    template <typename NodePtrIterator>  
    cutinfo_ptr runPartitionedSection(const NodePtrIterator& start, 
                                      const NodePtrIterator& end, uint key) {

      checkPartitionedSection(start, end, key);
      
      Base::enableChecks();

      Base::runSection(start, end, key);

      bool any_on = false;
      bool any_off = false;

      for(NodePtrIterator it = start; it != end; ++it) {
        if( (*it)->state() ) {
          any_on = true;
        } else {
          any_off = true;
        }
        if(any_on && any_off)
          break;
      }

      cutinfo_ptr cut = cutinfo_ptr(new CutInfo);

      if(any_on && !any_off) {
        for(NodePtrIterator it = start; it != end; ++it) {
          node_ptr n = *it;    
          n->template flipNode<1>(Base::lattice);
        }
        cut->any_cut = false;

#ifndef NDEBUG        
        for(NodePtrIterator it = start; it != end; ++it) {
          node_ptr n = *it;
          assert(n->matchesKey(key));
          assert(n->height == 0);
          assert(n->state() == 0);
        }
#endif

        return cut;

      } else if (any_off && !any_on) {
        cut->any_cut = false;


#ifndef NDEBUG        
        for(NodePtrIterator it = start; it != end; ++it) {
          node_ptr n = *it;
          assert(n->matchesKey(key));
          assert(n->height == 0);
          assert(n->state() == 0);
        }
#endif

        return cut;
      }
      
      cut->any_cut = true;
      
      cut->partitions[0] = partitioninfo_ptr(new PartitionInfo);
      cut->partitions[0]->is_on = false;

      cut->partitions[1] = partitioninfo_ptr(new PartitionInfo);
      cut->partitions[1]->is_on = true;

      cut->cut_value = 0;

      assert(cut->partitions[0]->nodes.empty());
      assert(cut->partitions[1]->nodes.empty());

      // First go through and set the state to the proper node.  all
      // these are currently eligible
      for(NodePtrIterator it = start; it != end; ++it) {
        node_ptr n = *it;    

        assert(n->matchesKey(key));

        bool n_on = n->state();
        
        cut->partitions[n_on ? 1 : 0]->nodes.push_back(n);

        uint key_state_check = Node::makeKeyState(key, 1);

        // Now see about the cut.
        if(!n_on) {
          for(uint ei = 0; ei < kernel_size; ++ei) {
            node_ptr nn = n + Base::step_array[ei];
          
            if(nn->_isKeyState(key_state_check)) {
              assert(nn->matchesKey(key) && nn->state());
              dtype cc = n->edgeCapacityFromSaturatedOffNodeToOnNode(ei);
              cut->cut_value += cc;
              cut->cut_edges.push_back(make_pair(nn, Base::reverseIndex(ei)));
            } else {
              assert(!(nn->matchesKey(key) && nn->state()));
            }
          }
        }

        // // Go through and fill in the values for the different partitions.
        // for(uint ei = 0; ei < kernel_size; ++ei) {
        //   node_cptr nn = n + Base::step_array[ei];          
          
        //   bool nn_on = nn->state();

        //   if(nn->matchesKey(key)) {
        //     if(!n_on && nn_on) {
        //       dtype cc = n->edgeCapacityFromSaturatedOffNodeToOnNode(ei);
        //       cut->cut_value += cc;
        //       cut->cut_edges.push_back(make_pair(n, ei));
        //     } 
        //   } else {
        //     assert(!nn->state());
            
        //     cut->gamma_sum[n_on] += n->influenceOnReduction(nn, ei);
        //     cut->qii_sum[n_on] += n->fv();
        //   }
        // }
      }

      for(node_ptr n : cut->partitions[0]->nodes) {
        assert(n->state() == 0);
      }
      
      // Flip all of these nodes back 
      for(node_ptr n : cut->partitions[1]->nodes) {
        n->template flipNode<1>(Base::lattice);
        assert(n->state() == 0);
      }

#ifndef NDEBUG        
        for(NodePtrIterator it = start; it != end; ++it) {
          node_ptr n = *it;
          assert(n->matchesKey(key));
          assert(n->height == 0);
          assert(n->state() == 0);
        }
#endif

      return cut;
    }

    void applyPartioningCut(cutinfo_ptr cut, uint key) {
      for(auto cutedge : cut->cut_edges) {
        node_ptr n = cutedge.first;
        uint ei = cutedge.second;
        node_ptr dest = n + Base::step_array[ei];
        assert(n->matchesKey(key));
        assert(dest->matchesKey(key));
        n->template saturate<0>(Base::lattice, dest, ei);
      }
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Method for walking the lattice and picking out sections.  This
    // can only be called from an external function; cannot be done while solving

    template <typename UnaryNodeCriterionFunction> 
    vector<node_ptr>  
    walkConnectedRegion(node_ptr seed_node, const UnaryNodeCriterionFunction& test_f) {
      
      if(DEBUG_MODE && ENABLE_EXPENSIVE_CHECKS) {
        for(const auto& n : Base::lattice) {
          assert_equal(n.height, 0);
        }
      }
      
      vector<node_ptr> nodes(1, seed_node);

      size_t visit_index = 0;

      do {
        node_ptr n = nodes[visit_index];

        for(uint ei = 0; ei < kernel_size; ++ei) {

          node_ptr nn = n + Base::step_array[ei];
          
          if(nn->height == 0
             && (Base::pushCapacity(n, nn, ei) 
                 + Base::pushCapacity(nn, n, Base::reverseIndex(ei)) > 0)
             && (test_f(nn))) {
           
            nn->height = 1;
            nodes.push_back(nn);
          }
        }
      } while( ( (++visit_index) != nodes.size()));

      for(node_ptr n : nodes) 
        n->height = 0;

      if(DEBUG_MODE && ENABLE_EXPENSIVE_CHECKS) {
        for(const auto& n : Base::lattice) {
          assert_equal(n.height, 0);
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
          
        if(Base::pushCapacity(n, nn, ei) + Base::pushCapacity(nn, n, Base::reverseIndex(ei)) > 0) 
          return false;
      }

      return true;
    }

    template <typename RegisterNodePair> 
    void constructNeighborhoodKeyPairs(RegisterNodePair& pair_reg) const {

      unordered_set<uint64_t> seen_keys;

      auto gen_key = [](uint key1, uint key2){
        return (uint64_t(min(key1, key2)) << 32) | (uint64_t(max(key1, key2)));
      };

      for(node_ptr n = Base::lattice.begin(); n != Base::lattice.end(); ++n) {
        seen_keys.insert(gen_key(n->key(), n->key()));

        for(uint ei = 0; ei < Base::lattice.kernel_positive_size; ++ei) {
          node_ptr nn = n + Base::step_array[ei];

          if(Base::pushCapacity(n, nn, ei) 
             + Base::pushCapacity(nn, n, Base::reverseIndex(ei)) > 0) {
            auto ret = seen_keys.insert(gen_key(n->key(), nn->key()));

            if(ret.second) { 
              // The element was inserted; hasn't been seen before
              pair_reg(n, nn);
            }
          }
        }
      }
    }

    template <typename ForwardIterator>
    set<uint> getNeighborhoodKeySet(const ForwardIterator& start, 
                                    const ForwardIterator& end, 
                                    uint this_key) const {
      set<uint> keys;

      for(ForwardIterator it = start; it != end; ++it) {
        node_ptr n = (*it);

        assert_equal(this_key, n->key());
        
        for(uint ei = 0; ei < Base::lattice.kernel_size; ++ei) {
          node_ptr nn = n + Base::step_array[ei];

          if(!nn->matchesKey(this_key) &&
             (Base::pushCapacity(n, nn, ei) 
              + Base::pushCapacity(nn, n, Base::reverseIndex(ei)) > 0)) {
            keys.insert(nn->key());
          }
        }
      }
      
      return keys;
    }

    template <typename ForwardIterator>
    void checkKeyIsGone(const ForwardIterator& start, 
                        const ForwardIterator& end, 
                        uint key) const {
      if(DEBUG_MODE) {
        for(ForwardIterator it = start; it != end; ++it) {
          node_ptr n = (*it);

          assert_notequal(key, n->key());
        }
      }
    }

    template <typename NodePtrIterator> 
    inline dtype getExcessInRegion(const NodePtrIterator& start, 
                                   const NodePtrIterator& end) {
      
      dtype flow_excess = 0;
      
      for(NodePtrIterator it = start; it != end; ++it) {
        node_ptr n = *it;
        flow_excess += Base::excess(n);
      }

      return flow_excess;
    }

    template <typename NodePtrIterator>  
    inline void setRegionToLambda_reference(const NodePtrIterator& start, 
                                            const NodePtrIterator& end, 
                                            dtype lambda) {
      
      // First, get the accurate offsets
      comp_type fv_avg = 0;
      comp_type qii_sum = 0; 
      size_t n_sum = 0;

      for(auto it = start; it != end; ++it) {
        node_ptr n = (*it);

        n->setOffsetAndScale(Base::lattice, 0, lambda); 

        fv_avg += n->r();
        ++n_sum;
      }

      dtype fv_offset = floorAverage(fv_avg, n_sum);

      for(auto it = start; it != end; ++it) 
        (*it)->setOffset(Base::lattice, fv_offset);
   
    }

    template <typename NodePtrIterator>  
    inline dtype setRegionToLambda(const NodePtrIterator& start, 
                                  const NodePtrIterator& end, 
                                  dtype lambda, bool pull_to_center = false) {

      comp_type zero_reference = 0;
      comp_type qii_sum = 0; 
      comp_type gamma_sum = 0; 
      
      long partition_size = 0;

      for(NodePtrIterator it = start; it != end; ++it) {
        node_ptr n = (*it);
        qii_sum += n->qii();
        gamma_sum += n->influence();
        ++partition_size;
      }

      // Now the base value should be the sum of all of these,
      // multiplied together.
      dtype fv_offset = floorAverage(Node::multFVScale(qii_sum, lambda) + gamma_sum, partition_size);

#ifndef NDEBUG
      {
        // First, get the accurate offsets
        comp_type fv_avg = 0;
        comp_type qii_sum = 0; 
        size_t n_sum = 0;

        for(auto it = start; it != end; ++it) {
          node_ptr n = (*it);

          n->setOffsetAndScale(Base::lattice, 0, lambda); 

          fv_avg += n->r();
          ++n_sum;
        }

        dtype _fv_offset = floorAverage(fv_avg, n_sum);
        assert_close(_fv_offset, fv_offset, 1);
      }
#endif

      for(NodePtrIterator it = start; it != end; ++it) {
        node_ptr n = (*it);
        // Go through and set the offset of each node such that it is
        // pulled towards the partition being nicer

        dtype offset = fv_offset + ((pull_to_center && (n->fv(lambda) > fv_offset)) ? 1 : 0);
        n->setOffsetAndScale(Base::lattice, offset, lambda);
        n->_debug_checkLevelsetMethodsNode(Base::lattice);
      }

      return fv_offset;
    }
  };
}

#include "../common/debug_flymake_test.hpp"

#endif /* _TV_PUSH_RELABEL_H_ */
