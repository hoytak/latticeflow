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
    };

    typedef shared_ptr<CutInfo> cutinfo_ptr;

    template <typename NodePtrIterator>  
    cutinfo_ptr runPartitionedSection(const NodePtrIterator& start, 
                                      const NodePtrIterator& end, uint key) {

#ifndef NDEBUG        
      for(NodePtrIterator it = start; it != end; ++it) {
        node_ptr n = *it;
        assert(n->matchesKey(key));
        assert(n->height == 0);
        assert(n->state() == 0);
      }
#endif

      Base::runSection(start, end, key);

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
      
      cut->partitions[0] = partitioninfo_ptr(new PartitionInfo);
      cut->partitions[0]->is_on = false;

      cut->partitions[1] = partitioninfo_ptr(new PartitionInfo);
      cut->partitions[1]->is_on = true;

      cut->cut_value = 0;

      // First go through and set the state to the proper node.  all
      // these are currently eligible
      for(NodePtrIterator it = start; it != end; ++it) {
        node_ptr n = *it;    

        assert(n->matchesKey(key));
        
        bool is_on = n->state();
        
        cut->partitions[is_on ? 0 : 1]->nodes.push_back(n);

        // Now see about the cut.
        for(uint ei = 0; ei < kernel_size; ++ei) {
          node_ptr nn = n + Base::step_array[ei];
          
          if(nn->matchesKey(key) && nn->state() != is_on) {
            cut->cut_value += capacityOfSaturated(n, nn, ei); 
            cut->cut_edges.push_back(make_pair(n, ei));
          }
        }
      }
      
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
      
      if(DEBUG_MODE) {
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

      if(DEBUG_MODE) {
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

    template <typename ForwardIterator, typename RegisterNode> 
    void constructNeighborhoodSet(const ForwardIterator& start, const ForwardIterator& end, 
                                  uint key, const RegisterNode& reg_node) const {

      for(ForwardIterator it = start; it != end; ++it) {
        node_ptr n = (*it);
        
        for(uint ei = 0; ei < Base::lattice.kernel_size; ++ei) {
          node_ptr nn = n + Base::step_array[ei];

          if(Base::pushCapacity(n, nn, ei) 
             + Base::pushCapacity(nn, n, Base::reverseIndex(ei)) > 0) {
            reg_node(nn);
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

        n->setOffsetAndScale(Base::lattice, 0, lambda); 

        fv_avg += n->cfv_predict();
        ++n_sum;
      }

      dtype fv_offset = ceilAverage(fv_avg, n_sum);

      for(auto it = start; it != end; ++it) {
        (*it)->setOffset(Base::lattice, fv_offset);
      }
    }
  };

}; 

#ifdef EMACS_FLYMAKE

#include "../kernels/kernels.hpp"
#include "../lattices/kernellattice.hpp"
#include "tv_flow_node.hpp" 

namespace latticeQBP {
  typedef KernelLattice<TVFlowNode<Star2d_4, long>, 2, Star2d_4> _TV_PRFlow_TestLattice;
  typedef TV_PRFlow<long, _TV_PRFlow_TestLattice> _TV_PRFlow_Test;

  template class TV_PRFlow<long, _TV_PRFlow_TestLattice>;
};

#endif



#endif /* _TV_PUSH_RELABEL_H_ */
