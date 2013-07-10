#ifndef _TV_PUSH_RELABEL_H_
#define _TV_PUSH_RELABEL_H_

#include "../common.hpp"
#include "../network_flow/push_relabel.hpp"

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

  protected:
    partitioninfo_ptr _addPartition(node_ptr n, uint key) {
      PartitionInfo* pi = new PartitionInfo;

      bool is_on              = pi->is_on     = n->on();
      dtype& cut_value        = pi->cut_value = 0;
      vector<node_ptr>& nodes = pi->nodes;

      nodes.push_back(n);

      for(size_t n_idx = 0; n_idx < nodes.size(); ++n_idx) {
        node_ptr n = nodes[n_idx];

        for(uint ei = 0; ei < kernel_size; ++ei) {
          node_ptr nn = n + Base::step_array[ei];
          if(nn->matchesKey(key)) {
            if(nn->state() == is_on) {
              nodes.push_back(nn);
            } else {
              if(DEBUG_MODE) {
                if(is_on) {
                  assert_lt(n->level(), nn->level());
                }
              }

              cut_value += capacityOfSaturated(n, nn, ei); 
            }
          }
        }

        n->clearKey();
      }

      return partitioninfo_ptr(pi);
    }

  public:
    template <typename NodePtrIterator>  
    vector<partitioninfo_ptr> cleanupAndGetBlocks(
       const NodePtrIterator& start, const NodePtrIterator& end, uint key = 0) {

      if(DEBUG_MODE) {
        for(NodePtrIterator it = start; it != end; ++it) {
          node_ptr n = *it;
          assert(n->matchesKey(key));
        }
      }

      vector<partitioninfo_ptr> piv;

      // First go through and set the state to the proper node.  all
      // these are currently eligible
      for(NodePtrIterator it = start; it != end; ++it) {
    
        // Set these nodes to the given key and partition
        node_ptr n = *it;
        if(!n->matchesKey(key))
          continue;
        else
          piv.push_back(_addPartition(n, key));
      }

      for(auto it : piv) {
        if(it->is_on != 0) {
          for(node_ptr n : it->nodes) {
            n->template flipNode<1>(Base::lattice);
          }
        }
      }
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

    // template <typename RegisterNodePair> {
      

    // }
  };

}; 

#endif /* _TV_PUSH_RELABEL_H_ */
