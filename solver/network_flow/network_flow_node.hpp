// Only meant to be included from within energy.hpp
#ifndef _NETWORK_FLOW_NODE_H_
#define _NETWORK_FLOW_NODE_H_

#include "../common.hpp"

namespace latticeQBP {

  ////////////////////////////////////////////////////////////////////////////////
  // BASE POLICY 

  struct NFNodeDefaultPolicy {
    static constexpr bool on_by_reduction = false;
    static constexpr bool enable_key_partitioning = false;
    static constexpr bool adjustment_mode = true;
  };

  ////////////////////////////////////////////////////////////////////////////////
  // NetworkFlowNode

  template <typename Kernel, typename dtype, typename Policy>
  class NetworkFlowNode {

  public:
    static constexpr bool reduction_mode  = Policy::on_by_reduction;
    static constexpr bool adjustment_mode = Policy::adjustment_mode;
    static constexpr bool enable_keys     = Policy::enable_key_partitioning;

    static constexpr bool simple_mode     = !(reduction_mode || adjustment_mode || enable_keys);

  public:
    NetworkFlowNode()
      : key_state(0)
      , height(0)
      , level_index(0)
      , reduction(0)
      , alpha(0)
#ifndef NDEBUG
      , q(0)
      , edges(0)
      , reduction_shift(0)
#endif
    {}

    static constexpr int n_edges = Kernel::size; 

    //////////////////////////////////////////////////
    // Helper methods that tie in other things...
    inline bool on() const { 
      if(reduction_mode)
        return reduction < 0;
      else
        return state() != 0;
    }

    inline unsigned int state() const {
      return enable_keys ? ((key_state & 0x1) != 0) : (key_state != 0);
    }

    template <int partition>
    static inline unsigned int makeKeyState(unsigned int key) {
      return ((key + 1) << 1) | partition;
    }

    template <int partition, typename Lattice> 
    inline unsigned int setKeyState(Lattice& lattice, unsigned int key) {
      static_assert(partition == 0 || partition == 1,
                    "Partition is the reference state; must be either 0 or 1.");

      if( (key_state & 0x1) != partition)
        this->template flipNode<1 - partition>(lattice);

      assert_equal((key_state & 0x1), partition);

      key_state = makeKeyState<partition>(key);
      
      return key_state;
    }
    
    inline bool _isKeyState(unsigned int _key_state) const {
      return key_state == _key_state;
    }

    template <int partition>
    inline void setKey(unsigned int key) {
      assert_equal(state(), partition);
        
      key_state = makeKeyState<partition>(key);
    }

    inline void clearKey() {
      key_state &= 0x1;
    }

    inline uint key() const {
      return (key_state >> 1) - 1;
    }

    inline bool keyIsClear() const {
      return ((key_state & (~uint(0x1))) == 0);
    }

    inline bool matchesKey(unsigned int key) const {
      return (key_state >> 1) == key + 1;
    }

  protected: 
    unsigned int key_state;

  public:

    // The height for the push-relabel thing
    typedef unsigned int level_index_type;

    level_index_type height;
    unsigned int level_index;

    // This also doubles as the capacity from source (positive) or cap to sink (negative)
    dtype reduction;

  protected:
    Array<dtype, n_edges> alpha;

#ifndef NDEBUG
    dtype q;
    Array<dtype, Kernel::size> edges;
  public:
    dtype reduction_shift;
  protected:
#endif

    template <int partition> inline void checkPartitioning() const {
      static_assert( adjustment_mode || partition == 0,
                     "Partition-based excess only allowed in non-simple mode.");

      if(adjustment_mode) 
        assert_equal(partition, key_state & 0x1);
    }

  public:
  
    ////////////////////////////////////////////////////////////////////////////////
    // Note that excess, in this case, is defined in relation to the
    // push-relabel algorithm; so excess-flow corresponds to deficit in
    // the reduction.  

    template <int partition=0> inline dtype excess() const {
      checkPartitioning<partition>();

      return ((partition == 1) ? 1 : -1) * (reduction);
    }

    template <int partition=0> inline dtype sinkFlow() const {
      checkPartitioning<partition>();
      return -excess<partition>();
    }

    template <int partition=0> inline dtype pushCapacity(int ei) const {
      checkPartitioning<partition>();

      if(adjustment_mode) {
        return (partition == 0) ? alpha[ei] : -alpha[ei];      
      } else {
        return alpha[ei];
      }
    }

    // This function only works if it is known the other node is in another region
    inline dtype capacityOfSaturated(int ei) const {
      dtype v = abs(alpha[ei]);

      if(DEBUG_MODE)
        assert_equal(v, edges[ei]);

      return v;
    }

    template <int partition=0> inline bool pushSaturated(int ei) const {
      return (pushCapacity<partition>(ei) <= 0);
    }

    void adjustReduction(dtype amount) {
      reduction += amount;
#ifndef NDEBUG
      reduction_shift += amount;
#endif
    }

    dtype level() const {
      assert_equal(reduction_shift, 0);
      return simple_mode ? reduction : reduction >> 1;
    }

    ////////////////////////////////////////////////////////////////////////////////
  protected:

    template <typename Lattice>
    void _debugVerifyNodeConsistency(Lattice& lattice, bool check_neighbors = true) { 

#ifndef NDEBUG

      if(!(this >= &(*lattice.begin())) && this < &(*lattice.end()))
        return;

      // Make sure the gain is what we'd like
      
      // dtype r_compare = 2*q;

      for(size_t i = 0; i < lattice.kernelSize(); ++i) {
        const auto& nn = lattice.neighbor(this, i);

        if(!lattice.isValidNode(nn))
          continue;

        uint rev_idx = lattice.reverseIndex(i);

        assert_equal(i, lattice.reverseIndex(rev_idx));

        assert(this == lattice.neighbor(nn, rev_idx));

        // r_compare += 2*edges[i] - abs(alpha[i]);

        // cout << "Testing nodes " << lattice._tag(this) << ", " << lattice._tag(nn) << ", ei = " << i << endl;

        if(simple_mode)
          assert_equal(abs(alpha[i]) + abs(nn->alpha[rev_idx]), edges[i]);
        else
          assert_equal(abs(alpha[i]) + abs(nn->alpha[rev_idx]), 2*edges[i]);

        if(simple_mode || adjustment_mode) {
          if(this->state() != 0 && nn->state() == 0) {
            assert_geq(alpha[i], 0);
            assert_leq(nn->alpha[rev_idx], 0);
          } 
          else if (this->state() == 0 && nn->state() != 0) {
            assert_leq(alpha[i], 0);
            assert_geq(nn->alpha[rev_idx], 0);
          }
        }
      }

      if(check_neighbors) {
        for(size_t i = 0; i < lattice.kernelSize(); ++i) {
          auto nn = lattice.neighbor(this, i);
          if(!lattice.isValidNode(nn))
            continue;
          nn->_debugVerifyNodeConsistency(lattice, false);
        }
      }

#endif
    }

  public:
  
    template <int start_state, typename LatticeType> inline void flipNode(LatticeType& lattice) {

#ifndef NDEBUG

      _debugVerifyNodeConsistency(lattice);

      for(size_t i = 0; i < lattice.kernel_size; ++i) 
        lattice.neighbor(this, i)->_debugVerifyNodeConsistency(lattice);
  
#endif

      // if(!simple_mode || DEBUG_MODE) {
      for(size_t i = 0; i < lattice.kernel_size; ++i)
        lattice.neighbor(typename LatticeType::value_direct_ptr(this), i)->alpha[lattice.reverseIndex(i)] *= -1;
      // }

      // All we need to do should be to set the lattice 
      assert(state() == 0 || state() == 1);

      if(enable_keys) {
        switch(start_state) {
        case 0: assert_equal(state(), start_state); key_state |= 1; break;
        case 1: assert_equal(state(), start_state); key_state &= (~((unsigned int)(1))); break;
        case 2: key_state ^= 0x1; break;
        default: break;
        }
      } else {
        switch(start_state) {
        case 0: assert_equal(state(), start_state); key_state = 1; break;
        case 1: assert_equal(state(), start_state); key_state = 0; break;
        case 2: key_state = !key_state; break;
        default: break;
        }
      }
#ifndef NDEBUG

      _debugVerifyNodeConsistency(lattice);

      for(size_t i = 0; i < lattice.kernelSize(); ++i) 
        lattice.neighbor(this, i)->_debugVerifyNodeConsistency(lattice);
#endif
    }

    template <int partition, typename Lattice, typename _node_ptr> 
    inline void pushExcess(Lattice& lattice, const _node_ptr& dest, unsigned int ei, dtype amount) {

#ifndef NDEBUG
      this->template checkPartitioning<partition>();
    
      dest->template checkPartitioning<partition>();

      assert(&(*dest) == &(*lattice.neighbor(this, ei)));
      assert_equal(dest->state(), partition);
      assert_leq(amount, pushCapacity<partition>(ei));
  
      _debugVerifyNodeConsistency(lattice);
      dest->_debugVerifyNodeConsistency(lattice);

      dtype src_pr_excess  = this->template excess<partition>();
      dtype dest_pr_excess = dest->template excess<partition>();

      for(size_t i = 0; i < lattice.kernel_size; ++i) {
        lattice.neighbor(this, i)->_debugVerifyNodeConsistency(lattice);
        lattice.neighbor(dest, i)->_debugVerifyNodeConsistency(lattice);
      }

#endif

      this->reduction -= ((partition == 1) ? 1 : -1) * amount;
      this->alpha[ei] -= amount;

      dest->reduction += ((partition == 1) ? 1 : -1) * amount;
      dest->alpha[Kernel::reverseIndex(ei)] += amount;
  
#ifndef NDEBUG
      _debugVerifyNodeConsistency(lattice);
      dest->_debugVerifyNodeConsistency(lattice);

      for(size_t i = 0; i < lattice.kernel_size; ++i) {
        lattice.neighbor(this, i)->_debugVerifyNodeConsistency(lattice);
        lattice.neighbor(dest, i)->_debugVerifyNodeConsistency(lattice);
      }

      dtype src_af_excess  = this->template excess<partition>();
      dtype dest_af_excess = dest->template excess<partition>();
  
      assert_equal(src_pr_excess - amount, src_af_excess);
      assert_equal(dest_pr_excess + amount, dest_af_excess);
#endif
    }

    template <int partition, typename Lattice, typename _node_ptr> 
    inline void saturate(Lattice& lattice, const _node_ptr& dest, unsigned int ei) {
      pushExcess<partition>(lattice, dest, ei, max(pushCapacity<partition>(ei), dtype(0)));
    }

    template <typename Lattice> 
    class NodeFiller {
    public:

      typedef typename Lattice::value_ptr node_ptr;

      NodeFiller(Lattice& _lattice) 
        : lattice(_lattice)
      {}

      template <typename NodePtr>
      inline void addE1(const NodePtr& n, dtype e0, dtype e1) const {
        dtype q = e1 - e0;

        assert(lattice.withinBounds(n));

        if(simple_mode)
          n->reduction += q;
        else 
          n->reduction += 2*q;

#ifndef NDEBUG
        n->_debugVerifyNodeConsistency(lattice);
        n->q += q;
#endif

      }

      template <typename NodePtr>
      inline void addE2(const NodePtr& n1, const NodePtr& n2, uint ei, 
                        dtype e00, dtype e01, dtype e10, dtype e11) const {

        // cout << "Adding nodes " << lattice._tag(n1) << ", " << lattice._tag(n2) << ", ei = " << ei << endl;

        dtype q1 = (e10 - e00);
        dtype q2 = (e01 - e00);
        dtype q12 = (e00 + e11 - e01 - e10);

        assert(lattice.withinBounds(n1));
        assert(lattice.withinBounds(n2));

        assert_leq(q12, 0);

        assert( (n1) <= (n2) );
      
        uint rev_idx = Kernel::reverseIndex(ei);

        assert(n2 == lattice.neighbor(n1, ei));
        assert(n1 == lattice.neighbor(n2, rev_idx));
        assert_equal(ei, lattice.reverseIndex(rev_idx));

        dtype e2 = -q12;

        assert_geq(e2, 0);

        if(simple_mode) {
          n1->reduction += q1 + q12;
          n2->reduction += q2;
          n1->alpha[ei] += e2;
        } else {
          n1->reduction += 2*q1 + q12;
          n2->reduction += 2*q2 + q12;
          n1->alpha[ei] += e2;
          n2->alpha[rev_idx] += e2;
        }

#ifndef NDEBUG
        n1->q += q1;
        n2->q += q2;

        n1->edges[ei] += e2;
        n2->edges[rev_idx] += e2;

        n1->_debugVerifyNodeConsistency(lattice);
        n2->_debugVerifyNodeConsistency(lattice);
#endif
      }
    private:
      Lattice& lattice;
    };
  };

};

#endif

