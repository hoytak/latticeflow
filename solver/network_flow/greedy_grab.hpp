
#ifndef _PRFLOW_PRIVATE_MEMBER_
#error "File greedy_grab.hpp cannot be included directly."
#else

template <int alloc_depth, int push_to_tips, int allow_nonzero_heights, int only_pull_needed> 
class GreedyGrab {
public:
  GreedyGrab(PRFlow& _prf)
    : prf(_prf), e(_prf.e), lattice(_prf.lattice)
  {}

private:
  PRFlow& prf;
  LatticeEnergy& e;
  Lattice& lattice; 

  Array<node_ptr, alloc_depth + 1> node_array;
  Array<DirDim, alloc_depth + 1> dd_array;
  Array<int, alloc_depth + 1> index_array;
  Array<ftype, alloc_depth +1> flow_request;

  unsigned int depth, max_depth, root_height;

  inline unsigned int requiredHeight(unsigned int _depth) const {
    if(!allow_nonzero_heights)
      return 0;
    else
      return (root_height > _depth) ? (root_height - _depth) : 0;
  }

  inline void tagNode(node_ptr n) const {
    if(!allow_nonzero_heights || n->height == 0) {
      assert_equal(n->height, 0);
      --(n->height);
    }
  }

  inline void untagNode(node_ptr n) const {
    if(!allow_nonzero_heights || n->height == (unsigned int)(-1))
      n->height = 0;
  }

  inline bool eligible(node_ptr n) const {
    return (n->on() == partition) && (n->height == requiredHeight(depth));
  }

  inline bool zeroHeight(int _depth) const {
    return (!allow_nonzero_heights || requiredHeight(_depth) == 0);
  }

#ifndef NDEBUG
  void checkConsistency() {
    // Make sure the flow requests and stuff are good

    for(unsigned int i = 1; i <= depth; ++i) {
      node_ptr root = node_array[i-1];
      node_ptr tip = node_array[i];
      DirDim dd = dd_array[i];

      assert(tip == lattice.neighbor(root, dd) );

      ftype cap = prf.pushCapacity(root, tip, dd);

      if(i != depth) {
	assert_gt(flow_request[i], 0);
	assert_gt(cap, 0);
      } else {
	assert_geq(flow_request[i], 0);
	assert_geq(cap, 0);
      }

      assert_leq(flow_request[i], cap);
    }
  }
#else
  inline void checkConsistency() {}
#endif

  // Advance index; returns true when things are completely done.
  template <int init_mode> inline bool _advanceIndex() {

  backup_depth:

    assert(depth >= 1);
    assert(depth <= max_depth);
      
    int& index = index_array[depth];

    if(init_mode) {
      index = 0;
    } else {
      untagNode(node_array[depth]);
      ++index;
    }

    assert_leq(index, lattice.neighborIteratorBound());

    node_ptr root = node_array[depth - 1];

    while(true) {

      if(index == lattice.neighborIteratorBound()) {
	if(init_mode) {
	  --depth;
	  return true;
	} else if (unlikely(--depth == 0) ) {
	  return true;
	} else {
	  goto backup_depth;
	}
      }

      DirDim dd = lattice.neighborDD(index);
      node_ptr n = lattice.neighbor(root, dd);

      ftype cap;

      if(eligible(n) && (cap = prf.pushCapacity(root, n, dd)) > 0) {

	tagNode(n);
	node_array[depth] = n;
	dd_array[depth] = dd;
	flow_request[depth] = min(cap, flow_request[depth-1] + prf.excess(root));
	return false;
      }

      ++index;
    }
        
    checkConsistency();
  }
  
  // Only returns true when we're completely done.
  inline bool advanceIndex() {
    bool r = _advanceIndex<0>();
    return r;
  }

  inline bool initNextLevel() {
    if(unlikely(depth == max_depth)) {
      return true;
    } else {
      ++depth;
      bool r = _advanceIndex<1>();
      return r;
    }
  }

  // Returns true if the edge was saturated 
  inline bool pullSRToRoot() {
    
    // Rules for pulling things: 

    // 1.  If flow_request[i] is not fulfilled, then the amount going
    // in is first allocated to flow_request[i-1] before what is
    // leftover is given to the local node.

    // 2. If the height of the node is greater than 0. then no new SF
    // can be created.  If the height is zero, then the flow is
    // governed by the only_pull_needed template parameter; if this is
    // true, then the flow is 

    if(push_to_tips) {
      assert_equal(depth, max_depth);
      assert_equal(prf.sinkFlow(node_array[depth]), 0);
    }

    ftype flow = push_to_tips ? flow_request[depth] : prf.sinkFlow(node_array[depth]);

    unsigned int saturated_tip = depth + 1;

    for(unsigned int i = depth; i != 0; --i) {
      // cout << "i = " << i << "; flow = " << flow << "; request = " << flow_request[i] 
      // 	   << "; cap = " << prf.pushCapacity(node_array[i-1], node_array[i], dd_array[i])
      // 	   << endl;

      if(i != 1)
	assert_leq(prf.sinkFlow(node_array[i-1]), 0);

      bool flow_sufficient = (flow >= flow_request[i]);
      if(flow_sufficient) {
	flow = flow_request[i];
	saturated_tip = i;
      }

      flow_request[i] -= flow;

      if(!only_pull_needed && zeroHeight(i-1)) {
	ftype cap = prf.pushCapacity(node_array[i-1], node_array[i], dd_array[i]);
	ftype desired = min(cap, prf.sinkFlow(node_array[i]));

	flow = max(flow, desired);
      }

      e.pushExcess<partition>(node_array[i-1], node_array[i], dd_array[i], flow);
    }

    bool has_saturated = (saturated_tip != depth + 1);

    if(has_saturated) {
      // Now pull back 
      for(unsigned int i = depth; i >= saturated_tip + 1; --i) {
	untagNode(node_array[i]);
      }
      depth = min(depth, saturated_tip);
    }

    checkConsistency();

    return has_saturated;
  }

public:
    
  // Returns true if the excess at the root has been satisfied
  bool run(node_ptr root, ftype requested_additional_sf_flow = 0, 
	   unsigned int _max_depth = alloc_depth) {

    assert(allow_nonzero_heights || root->height == 0 || root->height == 1);
    assert(!push_to_tips || root->height > _max_depth);
      
    max_depth = _max_depth;
    node_array[0] = root;
    tagNode(root);
    flow_request[0] = requested_additional_sf_flow;
    
    depth = 0;
    
    if(unlikely(initNextLevel())) {
      untagNode(root);
      return false;
    }

    while(true) {

      // INVARIANT: the tip is always a new one.
      node_ptr& tip = node_array[depth];

      // Can we send anything to the root? 
      if(!push_to_tips && prf.sinkFlow(tip) > 0) {
	bool edge_saturated = pullSRToRoot();

	if(unlikely(prf.sinkFlow(root) >= requested_additional_sf_flow)) {
	  goto flow_done;
	}
	  
	if(edge_saturated) {
	  if(unlikely(advanceIndex())) {
	    goto flow_done;
	  }
	  continue;
	}
      }

      // Now advance to the next one; returns true if there is no
      // way to advance

      if(initNextLevel()) {

	if(push_to_tips && depth == max_depth) {
	  bool edge_saturated = pullSRToRoot();

	  if(unlikely(prf.sinkFlow(root) >= requested_additional_sf_flow)) {
	    goto flow_done;
	  }

	  if(edge_saturated) {
	    if(unlikely(advanceIndex())) {
	      goto flow_done;
	    }
	    
	    continue;
	  }
	}

	// Okay, we're done
	if(unlikely(advanceIndex())) {
	  goto flow_done;
	}
      }
    }
      
  flow_done:
	
    for(unsigned int i = 0; i <= depth; ++i)
      untagNode(node_array[i]);
      
    bool result = (prf.sinkFlow(root) >= requested_additional_sf_flow);
    
    if(only_pull_needed && result) {
      assert_equal(prf.sinkFlow(root), requested_additional_sf_flow);
    }

    return result;
  }
};

#endif
