#ifndef _PUSH_RELABEL_H_
#define _PUSH_RELABEL_H_

#ifndef ENABLE_PR_CHECKS
#define ENABLE_PR_CHECKS false
#endif

#define NO_EXCESS    0
#define HAS_EXCESS   1
#define CHECK_EXCESS 2

#include <vector>
#include <set>
#include <queue>
#include <list>
#include <algorithm>
#include <iostream>
#include <memory>
#include <type_traits>

#include "../common.hpp"
#include "level_heap.hpp"
#include "fixed_sorting_functions.hpp"
#include "pr_auxilary_structures.hpp"
#include "../kernels/optimization_policies.hpp"

namespace latticeQBP {

  using namespace std;

  template <typename dtype, typename _KernelLattice, int partition> class PRFlow {
  protected:

    static_assert(partition == 0 || partition == 1, 
                  "partition must be 0 or 1.");

    static_assert(is_signed<dtype>::value, "dtype template parameter must be signed.");

    typedef _KernelLattice Lattice;
    typedef typename Lattice::Kernel Kernel;
    typedef KernelOptimizationPolicy<Kernel> OptPolicy; 
                  
    typedef typename Lattice::value_type Node;

    // Pull in the computation types 
    static constexpr bool reduction_mode  = Node::reduction_mode;
    static constexpr bool adjustment_mode = Node::adjustment_mode;
    static constexpr bool enable_keys     = Node::enable_keys;

  public:
    PRFlow(Lattice& _lattice)
      : lattice(_lattice)
      , key(0)
      , top_level(1)
      , step_array(lattice.jumpArray())
      , num_flips(0)
      , disable_printing(false)
      , enable_checks(false)
      , visited_nodes(OptPolicy::run_topology_restructure() 
                      ? 2*max(lattice.dimensions()) + 25 : 0)
      , level_heap(OptPolicy::run_topology_restructure() 
                   ? 2*max(lattice.dimensions()) + 25 : 0)
      , total_restructure_time(0)
      , node_run_counts(lattice.size(), 0)
    {
      levels.resize(2*max(lattice.dimensions()) + 25);

      if(DEBUG_MODE) {
        for(node_ptr n = lattice.begin(); n != lattice.end(); ++n) {
          n->setPosition(lattice.getCoords(n));
        }
      }
    }
 
  protected:
    typedef typename Lattice::value_ptr node_ptr;
    typedef typename Lattice::value_direct_ptr node_dptr;
    typedef typename Lattice::value_cptr node_cptr;
    typedef typename Lattice::value_type::level_index_type level_index_type;

    typedef unsigned int uint;
    static constexpr uint kernel_size = Lattice::kernel_size;
    static constexpr uint reverseIndex(uint idx) {return Lattice::reverseIndex(idx);}

    struct Level {
      vector<node_ptr> nodes;
      vector<node_ptr> active_nodes;
    };

    Lattice& lattice;
    vector<Level> levels;
    unsigned int key;

    size_t top_level;
    size_t top_level_with_excess;

    const Array<long, kernel_size> step_array;

    size_t num_flips;

    vector<node_ptr> node_buffer;

    bool disable_printing;

    bool enable_checks;

    ////////////////////////////////////////////////////////////////////////////////
    // Debug and consistency checks

  public:
    void disablePrinting() {
      disable_printing = true;
    }

  public:

    void _debug_forceVerifyAccurateLevel(size_t level, bool check_heights = false) {

      Level& lv = levels[level];

      set<node_ptr> all_nodes, seen_nodes;

      for(size_t i = 0; i < lv.nodes.size(); ++i) {
        node_ptr n = lv.nodes[i];

        assert(eligible(n));

        n->_debug_forceVerifyNodeConsistency(lattice, true);

        assert_equal(n->level_index, i);
        assert_equal(n->height, level);
        assert_equal(n->state(), partition);
        all_nodes.insert(n);

        if(check_heights) {
          for(uint j = 0; j < kernel_size; ++j) {
            node_ptr nn = n + step_array[j];

            if(eligible(nn) && pushCapacity(n, nn, j) > 0) {
              // if(!(nn->height >= level - 1))
              //  cout << "n = " << (n - lattice.begin()) << "nn = " << (nn - lattice.begin()) << endl;
              assert_geq(nn->height, level - 1);
            }
          }
        }
      }

      // for(size_t i = 0; i < lv.active_nodes.size(); ++i) {
      //   node_ptr n = lv.active_nodes[i];
      //   assert(seen_nodes.find(n) == seen_nodes.end());
      //   assert(all_nodes.find(n) != all_nodes.end());
      // }
    }

    void _debug_forceCheckPrerunStatus() {
      for(node_ptr n = lattice.begin(); n != lattice.end(); ++n) {
        n->_debug_forceVerifyNodeConsistency(lattice, true);

        if(!eligible(n))
          assert_equal(n->height, 0);

        assert_lt(n->height, levels.size());

        if(eligible(n) && excess(n) < 0) {
          if(n->height > 0)
            cerr << "  Examining Node  " << (n - lattice.begin()) << endl;
          assert_equal(n->height, 0);
        }
      }
    }

    void _debug_forceVerifyAll(bool check_all_nodes_in_level = true) {

      for(node_ptr n = lattice.begin(); n != lattice.end(); ++n) {
        n->_debug_forceVerifyNodeConsistency(lattice, true);

        if(!eligible(n))
          assert_equal(n->height, 0);

        assert_lt(n->height, levels.size());

        if(eligible(n) && n->height == 0)
          assert_leq(excess(n), 0);
      
        if(eligible(n) && excess(n) < 0) {
          if(n->height > 0)
            cerr << "  Examining Node  " << (n - lattice.begin()) << endl;
          assert_equal(n->height, 0);
        }
      }

      if(check_all_nodes_in_level) {
        for(size_t i = 1; i <= top_level; ++i) {

          Level& lv = levels[i];

          for(size_t j = 0; j < lv.nodes.size(); ++j) {
            node_ptr n = lv.nodes[j];

            n->height = 0;
          }
        }
   
        // Make sure all nodes are accounted for
        for(node_ptr n = lattice.begin(); n != lattice.end(); ++n)
          if(~(n->height))
            assert_equal(n->height, 0);
    
        for(size_t i = 1; i <= top_level; ++i) {

          Level& lv = levels[i];

          for(size_t j = 0; j < lv.nodes.size(); ++j) {
            node_ptr n = lv.nodes[j];

            n->height = i;
          }
        }
      }

      for(size_t i = 1; i <= top_level; ++i) {
        _debug_VerifyAccurateLevel(i, true);
        // assert_gt(levels[i].nodes.size(), 0);
      }
    }

    inline void enableChecks() { enable_checks = true; }
    inline void disableChecks() { enable_checks = false; }

#ifndef NDEBUG

    inline void _debug_checkPrerunStatus() { 
      if(ENABLE_PR_CHECKS || enable_checks) 
        _debug_forceCheckPrerunStatus();
    } 

    inline void _debug_VerifyAccurateLevel(size_t level, bool check = true) {
      if(ENABLE_PR_CHECKS || enable_checks) 
        _debug_forceVerifyAccurateLevel(level, check);
    }

    inline void _debug_VerifyAll(bool check = true) {
      if(ENABLE_PR_CHECKS || enable_checks) 
        _debug_forceVerifyAll(check);
    }
#else
    inline void _debug_checkPrerunStatus() { } 
    inline void _debug_VerifyAccurateLevel(size_t) { }
    inline void _debug_VerifyAll(bool = true) { }
#endif 

    inline void initStorage(size_t base_number) {

      size_t level_size = 4*base_number;
      size_t decrease = base_number;

      for(size_t i = 0; i != (levels.size() - 1) && levels[i].nodes.size() < level_size; ++i) {
        levels[i].nodes.reserve(level_size);
        levels[i].active_nodes.reserve(level_size);

        if(level_size <= decrease + levels[i + 1].nodes.size())
          break;
      
        level_size -= decrease;
      }
    }

    inline bool eligible(const node_ptr& n) const {
      return n->_isKeyState(key);
    }

    inline bool key_eligible(const node_ptr& n) const {
      return n->_isKey(key);
    }

    template <typename NodePtrIterator>
    void hotStart(const NodePtrIterator& start, const NodePtrIterator& end) {

      _debug_checkPrerunStatus();
      
#if HAVE_OUTPUT
      TimeTracker tt;
      tt.start();

      // // cerr << "Nodes: ";
      // for(NodePtrIterator it = start; it != end; ++it) {
      //   node_ptr n = lattice.resolve(*it);
      //   // cerr << lattice.index(n) << ",";
      // }
      // cerr << endl;
#endif 

      initQueue();

      // Hot start the nodes by building trees from all the sinks
      deque<node_ptr> original_nodes;

      for(NodePtrIterator it = start; it != end; ++it) {
      
        node_ptr n = lattice.resolve(*it);
        if(n->sinkFlow() > 0) {
          assert_equal(n->height, 0);
          original_nodes.push_back(n);
          n->height = 1;
        }
      }

      // Go through and add all these node neighbors to the queue
      deque<node_ptr> queue;

      for(node_ptr n : original_nodes) {
        for(uint i = 0; i < kernel_size; ++i) {

          node_ptr nn = n + step_array[i];

          if(eligible(nn) && nn->height == 0 
             && pushCapacity(nn, n, reverseIndex(i)) > 0) {
            addToLevel<CHECK_EXCESS, 0>(nn, 1);
            queue.push_back(nn);
          }
        }
      }
    
      // Now go through and run the queue
      size_t current_level = 0; // Triggers the resizing

      while(!queue.empty()) {
      
        node_ptr n = queue.front();
        assert(eligible(n));

        if(n->height != current_level) {
          ++current_level;
          assert_equal(current_level, n->height);
          assert_equal(queue.back()->height, n->height);

          if(unlikely(current_level >= levels.size())) {
            levels.resize(2*levels.size());
            cerr << "LEVELS RESIZE! " << endl;
          }

          levels[current_level].nodes.reserve(3*queue.size() / 2);
        }

        queue.pop_front();

        for(uint i = 0; i < kernel_size; ++i) {

          node_ptr nn = n + step_array[i];

          if(eligible(nn) && nn->height == 0 
             && pushCapacity(nn, n, reverseIndex(i)) > 0) {
            addToLevel<CHECK_EXCESS, 0>(nn, current_level + 1);
            queue.push_back(nn);
            assert_equal(nn->height, current_level + 1);
          }
        }
      }

      // Now, go through and see if the tree has missed any.  If so,
      // and these have excess, then they will end up flipping things later on
      deque<node_ptr> flip_nodes;

      for(NodePtrIterator it = start; it != end; ++it) {
      
        node_ptr n = lattice.resolve(*it);

        if(eligible(n) && n->height == 0 && excess(n) > 0) {
          flip_nodes.push_back(n);
          n->height = 1;
        }
      }
    
      if(unlikely(!flip_nodes.empty())) {
        while(!flip_nodes.empty()) {
          node_ptr n = flip_nodes.front();
          flip_nodes.pop_front();
          assert_equal(n->height, 1);
	
          for(uint i = 0; i < kernel_size; ++i) {
            node_ptr nn = n + step_array[i];

            if(pushCapacity(n, nn, i) > 0 && eligible(nn) && nn->height == 0) {
              flip_nodes.push_back(nn);
              nn->height = 1;
            }
          }

          flipNode(n);
        }
      }

      // Now, reset the original nodes to the right height
      for(node_ptr n : original_nodes) 
        n->height = 0;

#if HAVE_OUTPUT
      if(!disable_printing)
        cout << "Time taken in hotstart = " << tt.asString() << endl;
#endif
    }

    bool quickFlow() {
      // Go through and push stuff down hill if it's possible; kinda a
      // faster first pass with no theoretic value.

      for(size_t level = top_level; level != 1; --level) {
        _debug_VerifyAccurateLevel(level);
        _debug_VerifyAccurateLevel(level-1);
      
        Level& lv = levels[level];

        for(size_t node_index = 0; node_index < lv.active_nodes.size(); ++node_index) {
          node_ptr n = lv.active_nodes[node_index];

          if(!eligible(n) || n->height != level || excess(n) <= 0)
            continue;

          assert_gt(excess(n), 0);

          dtype remaining_excess = excess(n);

          for(uint i = 0; i < kernel_size; ++i) {
            node_ptr nn =  n + step_array[i];
            dtype cap;
            if( (cap = pushCapacity(n, nn, i)) != 0 && eligible(nn) && nn->height == level-1) {
              dtype amount = min(remaining_excess, cap);
              remaining_excess -= amount;
              push<1>(n, nn, i, amount);
	   
              if(remaining_excess == 0) {
                swap(lv.active_nodes[node_index], lv.active_nodes.back());
                lv.active_nodes.pop_back();
                --node_index;

                break;
              }
            }
          }
        }
      
        if(top_level_with_excess == level && lv.active_nodes.empty()) {
          --top_level_with_excess;
        }

        _debug_VerifyAccurateLevel(level);
        _debug_VerifyAccurateLevel(level-1);
      }

      ////////////////////////////////////////////////////////////
      // Now do the bottom level
      Level& lv = levels[1];

      for(size_t node_index = 0; node_index < lv.active_nodes.size(); ++node_index) {
        node_ptr n = lv.active_nodes[node_index];
        if(!eligible(n) || n->height != 1 || excess(n) <= 0)
          continue;

        dtype remaining_excess = excess(n);

        for(uint i = 0; i < kernel_size;++i) {
          node_ptr nn =  n + step_array[i];
          dtype cap, sink_flow;

          if( (cap = pushCapacity(n, nn, i)) != 0 
              && eligible(nn) 
              && nn->height == 0
              && (sink_flow = sinkFlow(nn)) > 0) {
	  
            dtype amount = min(min(remaining_excess, cap), sink_flow);
            remaining_excess -= amount;
            push<0>(n, nn, i, amount);
	  
            if(remaining_excess == 0) {
              swap(lv.active_nodes[node_index], lv.active_nodes.back());
              lv.active_nodes.pop_back();
              --node_index;
              break;
            }
          }
        }
      }

      if(top_level_with_excess == 1 && lv.active_nodes.empty()) {
        --top_level_with_excess;
        return true;
      }
    
      return false;
    }

    void clearLevels() {
      for(size_t level = 1; level <= top_level ; ++level) {

        _debug_VerifyAccurateLevel(level);

        for(size_t i = 0; i < levels[level].nodes.size(); ++i) {
          node_ptr n = levels[level].nodes[i];
          n->height = 0;
        }
      
        levels[level].nodes.clear();
        levels[level].active_nodes.clear();
      }

      top_level = 1;
      top_level_with_excess = 0;
    }

    //////////////////////////////////////////////////
    // The main queue structure

    inline bool nextInQueue(node_ptr& new_node, bool& base_level_empty) {

    restart:
      if(unlikely(top_level_with_excess == 0))
        return true;
 
      for(;levels[top_level_with_excess].active_nodes.size() == 0;) {
        --top_level_with_excess;
     
        if(unlikely(top_level_with_excess == 0))
          return true;
      }

      _debug_VerifyAccurateLevel(top_level_with_excess);
    
      Level& lv = levels[top_level_with_excess];

      // Remove the topmost node
      new_node = lv.active_nodes.back();
      lv.active_nodes.pop_back();
    
      // The second part is needed as the reshaping part doesn't remove
      // things from the active node queue.

      if( (OptPolicy::run_topology_restructure() 
           && (!eligible(new_node) 
               || new_node->height != top_level_with_excess))
          || excess(new_node) <= 0)
        goto restart;

      // if(new_node->level_index >= lv.nodes.size())  {
      //   cout << "MAJOR ERROR in LEVEL_INDEX -- out of range" << endl;
      //   cout << "level_index = " << new_node->level_index << "; size = " << lv.nodes.size() << endl;
      //   cout << "top_level_with_excess = " << top_level_with_excess << endl;
      //   cout << "active_nodes.size() = " << lv.active_nodes.size()  << endl;
      //   abort();
      // }

      // if(lv.nodes[new_node->level_index] != new_node){
      //   cout << "MAJOR ERROR in LEVEL_INDEX -- not same!!" << endl;
      //   abort();
      // }

      // Remove it from the level
      swap(lv.nodes[new_node->level_index], lv.nodes.back());
      lv.nodes[new_node->level_index]->level_index = new_node->level_index;
      lv.nodes.pop_back();
      base_level_empty = lv.nodes.empty();

      assert_gt(excess(new_node), 0);
      assert_equal(new_node->state(), partition);

      _debug_VerifyAccurateLevel(top_level_with_excess);
    
      return false;
    }

    template <int excess_status, int no_current_height_checks>
    inline void addToLevel(node_ptr n, size_t new_level) {
    
      const bool has_excess = ((excess_status == CHECK_EXCESS)
                               ? (excess(n) > 0)
                               : (excess_status == HAS_EXCESS));

      assert_geq(excess(n), 0); 

      if(excess_status == HAS_EXCESS)
        assert_gt(excess(n), 0);

      if(excess_status == NO_EXCESS)
        assert_equal(excess(n), 0);

      if(!no_current_height_checks) {
        assert_leq(n->height, top_level);
        assert(n->level_index >= levels[n->height].nodes.size() || 
               levels[n->height].nodes[n->level_index] != n);
      }

#ifndef NDEBUG
      unsigned int old_height = n->height;

      if(!no_current_height_checks)
        _debug_VerifyAccurateLevel(old_height);

      if(new_level <= top_level)
        _debug_VerifyAccurateLevel(new_level);
#endif

      if(new_level > top_level) {
        assert_equal(new_level, top_level + 1);
        ++top_level;
        if(unlikely(top_level >= levels.size())) {
          levels.resize(2*levels.size());
          cerr << "LEVELS RESIZE! " << endl;
        }
      }
    
      n->height = new_level;

      // Swap out the location of the 
      Level& lv = levels[new_level];
    
      size_t back_index = lv.nodes.size();
      lv.nodes.push_back(n);
      n->level_index = back_index;
    
      if(has_excess) {
        lv.active_nodes.push_back(n);

        if(top_level_with_excess < n->height)
          top_level_with_excess = n->height;
      }

#ifndef NDEBUG
      if(!no_current_height_checks)
        _debug_VerifyAccurateLevel(old_height);

      _debug_VerifyAccurateLevel(new_level);
#endif
    }

    inline void removeFromLevel(node_ptr n) {
      Level& lv = levels[n->height];

      // Remove it from the level
      swap(lv.nodes[n->level_index], lv.nodes.back());
      lv.nodes[n->level_index]->level_index = n->level_index;
      lv.nodes.pop_back();

      _debug_VerifyAccurateLevel(n->height);
    }

    inline void nodeNowHasExcess(node_ptr n) {
      assert_equal(n->state(), partition);
      assert_gt(excess(n), 0);
    
      levels[n->height].active_nodes.push_back(n);
      if(top_level_with_excess < n->height)
        top_level_with_excess = n->height;

      _debug_VerifyAccurateLevel(n->height);
    }

    inline void flipNode(node_ptr n, bool check_top_level = false) {
      if(check_top_level && levels[n->height].nodes.empty()) {
        if(unlikely(top_level > n->height)) {
          raiseAndFlipLevelAndAbove(n->height);
        } else {
          top_level = n->height - 1;
        }
      }

      assert_geq(excess(n), 0);
      n->height = 0;
      n->template flipNode<partition>(lattice);
      ++num_flips;
    }

    void raiseAndFlipLevelAndAbove(size_t level) {

      if(DEBUG_MODE) {
        for(size_t i = top_level + 1; i < levels.size(); ++i) {
          assert(levels[i].nodes.empty());
        }
      }
    
      // Pop off the queue
    
      for(size_t lv = level; lv <= top_level; ++lv) {
        for(node_ptr n : levels[lv].nodes) {
          flipNode(n);
        }
      
        levels[lv].nodes.clear();
      }

      top_level = level - 1;
    }

    void initNodeHeight(node_ptr n) {
      addToLevel<CHECK_EXCESS, 0>(n, 1);
    }

    void finish() {

      _debug_VerifyAll();

      assert(levels[0].nodes.empty());

      for(size_t level = 1; level <= top_level ; ++level) {

        for(size_t i = 0; i < levels[level].nodes.size(); ++i) {
          node_ptr n = levels[level].nodes[i];

          assert_leq(excess(n), 0);
          n->height = 0;
        }
      
        levels[level].nodes.clear();
      }

      for(size_t i = top_level + 1; i < levels.size(); ++i) {
        assert(levels[i].nodes.empty());
      }

      top_level = 1;
    }
  
    ////////////////////////////////////////////////////////////////////////////////
    // Debug stuff


    //////////////////////////////////////////////////////////////////////
    // Low level operAtions

    inline dtype pushCapacity(const node_cptr& src, const node_cptr& dest, uint ei) const {
      return max(src->template pushCapacity<partition>(ei), dtype(0));
    }

    inline dtype capacityOfSaturated(const node_cptr& src, const node_cptr& dest, uint ei) const {
      assert_notequal(src->state(), dest->state());
      src->_debugVerifyNodeConsistency(lattice);
      dest->_debugVerifyNodeConsistency(lattice);

      dtype ret = src->capacityOfSaturated(ei);

      assert_equal(ret, pushCapacity(dest, src, reverseIndex(ei)) + pushCapacity(src, dest, ei));

      return ret;
    }

    template <int node_may_have_excess>
    inline void push(const node_ptr& src, const node_ptr& dest, unsigned int ei, dtype amount) {

      if(unlikely(amount == 0))
        return;

      src->template pushExcess<partition>(lattice, dest, ei, amount);
    
      dtype dest_excess = excess(dest);

      if(node_may_have_excess && likely(dest_excess > 0 && dest_excess <= amount)) {
        nodeNowHasExcess(dest);
      }

      if(OptPolicy::run_topology_restructure() ){
        if(dest_excess >= 0 && amount > dest_excess) {
          node_run_counts[dest - lattice.begin()] 
            += OptPolicy::topology_restructure_sink_fill_bonus();
        }
      }

      if(!node_may_have_excess)
        assert_leq(excess(dest), 0);
    }

    inline dtype excess(const node_cptr& n) const {
      return n->template excess<partition>();
    }

    inline dtype sinkFlow(const node_cptr& n) const {
      return max(n->template sinkFlow<partition>(), dtype(0));
    }

    inline void addToQueue(node_ptr n) {
      if(likely(excess(n) > 0))
        addToLevel<HAS_EXCESS, 0>(n, 1);
    }


    ////////////////////////////////////////////////////////////////////////////////
    // Stuff to deal with the triggering for the node count

    vector<node_ptr> visited_nodes;
    PRLevelHeap<node_ptr> level_heap;
    TimeTracker total_restructure_time;
    vector<uint> node_run_counts;

    void addNodeRunCountEvent(node_ptr n) {
#ifndef NDEBUG
      size_t c = (++ node_run_counts[n - lattice.begin()]);
      assert(c <= OptPolicy::topology_restructure_trigger_threshold());
#else
      (++ node_run_counts[n - lattice.begin()]);
#endif
    }

    inline void clearNodeRunCounter(node_ptr n) {
      node_run_counts[n - lattice.begin()] = 0;
    }

    inline bool topologyRestructureTrigger(node_ptr n) const {
      return (node_run_counts[n - lattice.begin()] 
              >= OptPolicy::topology_restructure_trigger_threshold());
    }

    ////////////////////////////////////////////////////////////////////////////////
    //

    struct ExpansionNode {
      node_ptr node;
      level_index_type height;
      bool has_excess;
    };

    size_t restructureFromNode(node_ptr seed_node) {
    
      total_restructure_time.start();

      // Add everything to the level queue, flagging as active (n->height = max)
      level_heap.reset();
      visited_nodes.clear();

      _debug_VerifyAll(false);

      const level_index_type active_height = ~level_index_type(0);
      static_assert(~level_index_type(0) > level_index_type(0), 
                    "level_index_type not unsigned.");

      // Set up with the seed node.  Note that this has already been
      // removed from the level
      level_heap.push(seed_node->height, seed_node);
      seed_node->height = active_height;

      visited_nodes.push_back(seed_node);
      assert_gt(excess(seed_node), 0);

      size_t n_active_nodes_with_excess = (excess(seed_node) > 0) ? 1 : 0;

      // cerr << "Running restructuring from node " << (seed_node - lattice.begin()) << endl;

      // Time to run it
      level_heap.initIteration();

      Array<ExpansionNode, kernel_size> nl;

      while(!level_heap.empty()) { //  && n_active_nodes_with_excess != 0) {
      
        size_t level = level_heap.currentLevel();
        node_ptr n = level_heap.current();

        if(n->height != active_height) {

          // cerr << "Skipping visited node " << (n - lattice.begin()) << endl;
          if(excess(n) > 0)
            --n_active_nodes_with_excess;

          level_heap.pop();
          continue;
        }

        // cerr << "Running on node " << (n - lattice.begin()) << endl;
        // cerr << "   >>  height = " << level << endl;

        // Algorithm: Pick node off that are at the lowest level.  Add
        // all neighbors that are not downhill from it and that can be
        // pushed to.

        // We end things by not adding any new nodes to the queue and
        // running it all out when we don't have any active nodes with
        // excess.  

        level_index_type replace_height = active_height;
        uint nl_size=0;

        for(uint i = 0; i < kernel_size; ++i) {

          node_ptr nn = n + step_array[i];
	
          if(pushCapacity(n, nn, i) != 0 
             && eligible(nn)
             && nn->height != active_height) {

            nl[nl_size++] = {nn, nn->height, excess(nn) > 0};
            replace_height = min(replace_height, nn->height);
          }
        }

        // cout << "  replace_height = " << replace_height << endl;

        if(replace_height < level) {

          assert_leq(level, replace_height + 1);

          addToLevel<CHECK_EXCESS, 1>(n, replace_height + 1);

          // cerr << "Node " << (n - lattice.begin()) << " can push to neighbor; removing from active set. " << endl;

          if(excess(n) > 0) 
            --n_active_nodes_with_excess;

          // Add in all the nodes around this that are active so the
          // next rounds also build off this
          for(uint i = 0; i < kernel_size; ++i) {

            node_ptr nn = n + step_array[i];
	
            if(eligible(nn) 
               && nn->height == active_height 
               && pushCapacity(nn, n, reverseIndex(i)) > 0) {
	  
              // cerr << "Readding node " << (n - lattice.begin()) << " into active set. " << endl;

              level_heap.push(n->height + 1, nn);
	    
              // This must be done at all times, as 
              if(excess(nn) > 0) 
                ++n_active_nodes_with_excess;
            }
          }

        } else if(nl_size) {

          for(uint i = 0; i < kernel_size; ++i) {

            node_ptr nn = nl[i].node;
            level_index_type h = nl[i].height;
            bool has_excess = nl[i].has_excess;

            // This shouldn't be a sink
            assert_geq(excess(nn), 0);

            // cerr << "Adding new node " << (n - lattice.begin()) << " into active set. " << endl;
            removeFromLevel(nn);
            nn->height = active_height;
            level_heap.push(h, nn);
            visited_nodes.push_back(nn);

            if(has_excess)
              ++n_active_nodes_with_excess;

            if(!(--nl_size))
              break;
          }
        }

        // Go to the next one; have to advance at the end of the loop to
        // ensure things are not added to the current level (messing up internal stuff).
        level_heap.pop();
      }

      // Check to see if there is no sink to push the given node excess
      // to; in this case, the active set has just exausted itself.  
    
      // Note: is it possible that such a region is gauranteed to be
      // contiguous?  This way of doing things will have problems
      // otherwise...
      if(n_active_nodes_with_excess != 0) {
        assert(level_heap.empty());
      
        for(node_ptr n : visited_nodes) {
          if(n->height == active_height) {
            assert_geq(excess(n), 0);
            flipNode(n);
          }
          clearNodeRunCounter(n);
        }

      } else {

        while(!level_heap.empty()) {
          
          // size_t level = level_heap.currentLevel();
          node_ptr n = level_heap.current();

          if(n->height != active_height) {
            level_heap.pop();
            continue;
          }

          level_index_type replace_height = active_height;
      
          for(uint i = 0; i < kernel_size; ++i) {

            node_ptr nn = n + step_array[i];
	
            if(pushCapacity(n, nn, i) != 0 
               && eligible(nn)
               && nn->height != active_height) {

              replace_height = min(replace_height, nn->height);
            }
          }

          assert(replace_height != active_height);

          addToLevel<CHECK_EXCESS, 1>(n, replace_height + 1);

          // Add in all the nodes around this that are active so the
          // next rounds also build off this
          for(uint i = 0; i < kernel_size; ++i) {

            node_ptr nn = n + step_array[i];
	  
            if(eligible(nn) 
               && nn->height == active_height 
               && pushCapacity(nn, n, reverseIndex(i)) > 0 ) {
	  
              level_heap.push(replace_height + 1, nn);
            }
          }
        }
      }

      // Now, make sure there are none left that we haven't covered...
      for(node_ptr n : visited_nodes) {
        assert(n->height != active_height);
        clearNodeRunCounter(n);
      }

      _debug_VerifyAll();

      total_restructure_time.stop();

      return visited_nodes.size();
    }

    ////////////////////////////////////////////////////////////////////////////////

    static constexpr int kernel_compact_mode = OptPolicy::use_tiny_kernel_mode();
  
    typedef PossibleNeighbor<level_index_type, kernel_size, kernel_compact_mode> PNeighbor;

    void _run() {

      if(levels[1].nodes.empty())
        return;

      // Go through, popping off the highest one in the queue and
      // blasting through things until all is done!

      node_ptr n = lattice.end();

      Array<PNeighbor, kernel_size> hl;
      Array<node_ptr, kernel_size> nl;

      while(true) {

      next_iteration:

        bool base_level_empty = false;
        _debug_VerifyAll();
        bool done = nextInQueue(n, base_level_empty);
      
        if(unlikely(done))
          break;

        unsigned int base_height = n->height;
        dtype remaining_excess = excess(n);

        assert_gt(remaining_excess, 0);

        unsigned int first_nn = kernel_size;
        unsigned int last_valid = 0;

        for(uint i = 0; i < kernel_size; ++i) {

          node_ptr nn = nl[i] = n + step_array[i];

          dtype cap;
          bool valid = ( ( (cap = pushCapacity(n, nn, i)) != 0) && eligible(nn));
          //bool valid = (eligible(nn) && ((cap = pushCapacity(n, nn, i)) != 0));

          if(! OptPolicy::use_truncated_sorting()) {
            hl[i].set(valid ? (nn->height + 1) : 0, i);
	
            if(valid) {
              --first_nn;
              assert_leq(base_height, hl[i].height() );
              last_valid = i;
            }
          } else {
            if(valid) {
              hl[--first_nn].set(nn->height + 1, i);
              last_valid = i;
            }
          }
        }

        if(unlikely(first_nn == kernel_size)) {
          flipNode(n, base_level_empty);
          _debug_VerifyAll();
          goto next_iteration;
        }

        _debug_VerifyAccurateLevel(top_level_with_excess);

        if(first_nn == kernel_size - 1) {
          hl[first_nn] = hl[last_valid];
        } else {
          if(OptPolicy::use_fixed_sized_sorter()) {
            fixed_sorter<PNeighbor, kernel_size>::sort(hl.begin(), hl.end());
          } else {
            if(OptPolicy::use_truncated_sorting()) {
              sort(hl.begin() + first_nn, hl.end());
            } else {
              sort(hl.begin(), hl.end());
            }
          }
        }
#ifndef NDEBUG
        for(unsigned int i = first_nn; i < kernel_size; ++i)
          assert_geq(hl[i].height(), 1);
#endif

        if(base_height == 1 && hl[first_nn].nonzeroHeightEqualsOne()) {

          // Calculate if we can just send flow down there...
          for(unsigned int i = first_nn; 
              i < kernel_size && hl[i].nonzeroHeightEqualsOne(); ++i) {

            int index = hl[i].index();
            node_ptr nn = nl[index];
            dtype push_capacity = pushCapacity(n,nn,index);
            dtype amount = min(min(push_capacity, remaining_excess), sinkFlow(nn));

            remaining_excess -= amount;

            push<0>(n, nn, index, amount);

            if(push_capacity == amount) {
              // We saturated that edge, no longer need it for consideration
              swap(hl[i], hl[first_nn]);
              ++first_nn;
            }

            if(unlikely(remaining_excess == 0)) {
              addToLevel<NO_EXCESS, 0>(n, 1);

              _debug_VerifyAll();
              goto next_iteration;
            }
          }

          // Okay, now we need to init this level
          for(uint i = first_nn; i < kernel_size && hl[i].nonzeroHeightEqualsOne(); ++i) {
            assert_gt(pushCapacity(n, nl[hl[i].index()], hl[i].index()), 0);

            initNodeHeight(nl[hl[i].index()]);
            hl[i].incrementHeight();
          }
	
          base_level_empty = false;
        }
      
        for(unsigned int i = first_nn; i < kernel_size; ++i) {

          if(unlikely(base_level_empty && hl[i].heightIsGreaterThan(base_height))) {

            raiseAndFlipLevelAndAbove(base_height);

            // Now, can't push anywhere, so this node is done...
            flipNode(n);

            _debug_VerifyAll();
            goto next_iteration;
          }

          if(unlikely(OptPolicy::run_topology_restructure()) 
             && hl[i].heightIsGreaterThan(base_height)) {
            if(topologyRestructureTrigger(n)) {

              size_t n_hit = restructureFromNode(n);

              if(100*n_hit >= 
                 OptPolicy::quickflow_percent_after_topology_restructure_threshold() * lattice.size())

                quickFlow();

              _debug_VerifyAll();
              goto next_iteration;

            } else {
              addNodeRunCountEvent(n);
            }
          }

          int index = hl[i].index();
          const node_ptr& nn = nl[index];
          dtype amount = min(pushCapacity(n,nn,index), remaining_excess);

          remaining_excess -= amount;

          push<1>(n, nn, index, amount);

          if(remaining_excess == 0) {
            addToLevel<NO_EXCESS, 0>(n, hl[i].height() );

            assert(!levels[hl[i].height()].nodes.empty());

            _debug_VerifyAll();
            goto next_iteration;
          }
        }

        // Now, can't push anywhere, so this node is done...
        flipNode(n, base_level_empty);
        _debug_VerifyAll();
      }
      
      _debug_VerifyAll();
 
      // Now get rid of all the nodes on the queue
      finish();
    
      top_level = 1;
      top_level_with_excess = 0;
    
      // And we're done!
#ifndef NDEBUG
      bool fail = false;
      for(node_ptr n = lattice.begin(); n != lattice.end(); ++n) {
        if(n->height > 0) {
          cerr << "ERROR: Node " << (n - lattice.begin()) 
               << " does not have the height set right!!!" << endl;
          fail = true;
        }
      }
    
      assert(!fail);
#endif
    }

    ////////////////////////////////////////////////////////////////////////////////
    // higher level control structures

    void initQueue(size_t reserve_size = 0) {
      if(reserve_size != 0)
        initStorage(reserve_size);
      
      top_level = 1;
      top_level_with_excess = 0;
      
      num_flips = 0;
    }
  
  public:

    bool run() {
      _debug_checkPrerunStatus();

      // First do a simple iteration
      size_t starting_size = 0;
      size_t easy_count = 0;

      if(false && OptPolicy::init_hotstart()) {
        hotStart(lattice.begin(), lattice.end());
        if(OptPolicy::init_quickflow()) quickFlow();
        if(OptPolicy::init_rerun_hotstart()) {
          clearLevels();       
          hotStart(lattice.begin(), lattice.end());
        }
      } else {
        vector<node_ptr> active_nodes;

        for(node_ptr n = lattice.begin(); n != lattice.end(); ++n) {
          if(n->state() == partition && n->excess() > 0) {
            ++starting_size;
            active_nodes.push_back(n);
          }
        }
    
        initQueue(active_nodes.size());

        for(auto& n : active_nodes) {
          addToQueue(n);
        }
      }

      _debug_VerifyAll();

      _run();

      if(HAVE_OUTPUT && !disable_printing)
        cout << "Partition " << partition
             << ": Start size = " << starting_size << ", simple runner eliminated "
             << easy_count << "; " << num_flips << "/" << lattice.sizeWithinBounds()
             << " total flips. \n"
             << "Time spent in topology_restructure: " << total_restructure_time.asString() << "." << endl;

      _debug_VerifyAll();

      if(starting_size == 0)
        return true;
    
      return false;
    }

    ////////////////////////////////////////////////////////////////////////////////
    // stuff for working on specific sections.


    template <typename NodePtrIterator>  
    inline void prepareSection(const NodePtrIterator& start, 
                               const NodePtrIterator& end, 
                               uint _key = 0, bool do_state_cleaning = false) {

      // First go through and set the state to the proper node.  all
      // these are currently eligible
      for(NodePtrIterator it = start; it != end; ++it) {
    
        // Set these nodes to the given key and partition
        node_ptr n = *it;
        if(do_state_cleaning) {
          n->template setKeyState<partition>(lattice, _key);
        } else {
          assert_equal(n->state(), partition);
          n->template setKey<partition>(_key);
        }
      }
      // All that's needed here; just prepares things for runSection
    }

    template <typename NodePtrIterator>
    bool runSection(const NodePtrIterator& start, const NodePtrIterator& end, uint _key = 0) {
      // First do a simple iteration
      size_t set_size = 0;
      size_t starting_size = 0;

      key = Node::template makeKeyState<partition>(_key);

#ifndef NDEBUG
      {
        set<node_ptr> _nodes; 
        for(auto it = start; it != end; ++it)
          _nodes.insert(*it);
        
        for(node_ptr n = lattice.begin(); n != lattice.end(); ++n) {
          if(_nodes.find(n) == _nodes.end()) 
            assert(!eligible(n));
          else
            assert(eligible(n));
        }

        for(auto it = start; it != end; ++it) {
          node_ptr n = *it;
          assert_equal(n->state(), 0);
          assert_equal(n->key(), _key);
          assert(eligible(n));
        }
      } 
#endif 

      if(OptPolicy::init_hotstart()) {
        hotStart(start, end);
        if(OptPolicy::init_quickflow()) quickFlow();

        if(OptPolicy::init_rerun_hotstart()) {
          clearLevels();       
          hotStart(start, end);
        }
      } else {
        vector<node_ptr> active_nodes;

        for(NodePtrIterator it = start; it != end; ++it) {
          node_ptr n = *it;

          ++set_size;

          assert(eligible(n));

          if(n->excess() > 0) {
            ++starting_size;
            active_nodes.push_back(n);
          }
        }
    
        initQueue(active_nodes.size());

        for(auto& n : active_nodes) {
          addToQueue(n);
        }
      }

      _run();

#ifndef NDEBUG
      if(!disable_printing) {
        cerr << "Partition " << partition
             << ": Set size = " << set_size << "; queue start size = " 
             << starting_size << "; " << num_flips << " total flips. " << endl;
      }
#endif

      if(starting_size == 0)
        return true;
    
      return false;
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Additional utilities for some particular operations 

    template <typename NodePtrIterator>  
    inline void cleanupSection(const NodePtrIterator& start, 
                               const NodePtrIterator& end, 
                               uint _key = 0, bool do_state_cleaning = true) {

      // First go through and set the state to the proper node.  all
      // these are currently eligible
      for(NodePtrIterator it = start; it != end; ++it) {
    
        // Set these nodes to the given key and partition
        node_ptr n = *it;
        assert(n->matchesKey(_key));

        if(do_state_cleaning && n->state() != partition) 
          n->template flipNode<(1 - partition)>(lattice);

        assert_equal(n->state(), partition);
        n->clearKey();
      }
    }
  };
};
#endif
