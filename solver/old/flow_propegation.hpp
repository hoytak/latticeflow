// File meant to be included only from the flow propegation stuff

#ifndef _ENERGY_PRIVATE_MEMBER_
#error "File energy_debug.hpp cannot be included directly."
#else

/*********************************************************************************

Algorithm:

 

Each time a new level is added:

- See if the current excess is possibly sufficient to cover the need
  at the root node.  If so, then send it down.  However, propegate the
  mass back out if external need exists, and recallibrate the need at
  the center.  

- 

Notes:

- If you pull mass along a non-tree route, you have to adjust the
  local mass vs. the node excess(), any subtracted from one must be
  subtracted from the other.

- Pulling mass from a dup-node is allowed, up the mass at the node.
  This means that the need() at a node may be negative, with the
  available mass compensating.  In this case, the back-propegation
  needs to ensure that the need() is satisfied with the mass.

 ********************************************************************************/

struct DupNode;
typedef list<DupNode> dn_container;
typedef typename dn_container::iterator dn_ptr;

struct DupNode {
  node_ptr node;
  fn_ptr parent;
  DirDim step_to_parent;
  dn_ptr next_back_dup_node;
};

struct FlowNode {
  FlowNode() 
    : step_to_parent(DirDim(0,0))
    , depth(0)
    , first_parent_set(false)
    , inconsistent(false)
    , flow_to_parent(0)
    , need_upper_bound(0)
    , mass(0)
    , need_at_node_and_above(0)
    , external_need_at_node_and_above(0)
    , need_at_branch_tips(0)
  {}

  node_ptr node;
  fn_ptr parent, first_child;
  dn_ptr first_forward_dup_node;
  dn_ptr first_back_dup_node;
  DirDim step_to_parent;
  unsigned int depth;
  bool inconsistent;
  ftype flow_to_parent, need_upper_bound, mass;
  ftype need_at_node_and_above, external_need_at_node, external_need_at_node_and_above;
};

////////////////////////////////////////////////////////////////////////////////
// The container holding the tree; for easy traversal of things there

template <int partition> class Flow {
private:
  Lattice& lattice;

  fn_container container;
  fn_ptr next_in_queue;
  fn_ptr first_non_needy_neighborhood_node;

  dn_container dup_container;

  deque<fn_ptr> tree_tips;

  ftype total_unmet_need;

public:
  FlowNodeContainer(Lattice& _lattice, node_ptr root) 
    : lattice(_lattice)
    , next_in_queue(container.begin())
    , first_non_needy_neighborhood_node(container.end())
    , total_unmet_need(root.need())
  {
    container.push_back(FlowNode());
    container.back().node = root;
    container.back().parent = container.end();
    container.back().first_child = container.end();
    container.back().depth = 0;
    container.back().need_upper_bound = root->need();

    assert(total_unmet_need > 0);
  }

private:

  ////////////////////////////////////////////////////////////
  // Node query operations

  fn_ptr newFlowNode(node_ptr node, fn_ptr parent, const DirDim& step_to_child) {
    
    assert(!node->isVisited());
    assert(node->need() <= 0 || parent->node->need() > 0);

    fn_ptr r = container.insert((node->need() > 0)
				? first_non_needy_neighborhood_node 
				: container.end(),
				FlowNode());

    r->node			= node;
    r->parent			= parent;
    r->first_child		= container.end();
    r->first_forward_dup_node   = dup_container.end();
    r->first_back_dup_node      = dup_container.end();
    r->depth			= parent->depth + 1;
    r->step_to_parent		= r->step_to_child.reversed();
    r->need_upper_bound		= (min(parent->need_upper_bound, 
				       flowCapacity(node, _parent->node, _step_to_parent))
				   + node->need());
    r->mass			= node->excess();
    
    assert(r->need_upper_bound > 0);
    
    if(node->need() > 0)
      assert_gt(parent->node->need(), 0);

    node->setVisitor(r);

    return r;
  }
  
  inline dn_ptr newDupNode(node_ptr node, fn_ptr parent, const DirDim& step_to_child) {
    dn_ptr dn = dup_container.insert(dup_container.end(), DupNode());

    dn->node = node;
    dn->parent = parent;
    dn->step_to_parent = step_to_child.reversed(); 

    dn->next_back_dup_node = node->visitor()->dup_node_chain;

    node->visitor()->dup_node_chain = dn;
    
    return dn;
  }

  inline fn_ptr popQueue() {
    return (next_in_queue++);
  }

  ////////////////////////////////////////////////////////////
  // Querying operations

  bool nodeIsRoot(fn_ptr fn) const {
    return (fn->depth == 0);
  }

  bool nodeIsConsistent(fn_ptr fn) const {
    if(fn->consistent) {
      assert_equal(fn->excess, fn->node->excess());
      assert_equal(fn->flow_to_parent, 0);
    }
    
    return fn->consistent;
  }

  ////////////////////////////////////////////////////////////
  // movement operations 

  void sendMassTowardsRoot(fn_ptr fn) {
    assert(! nodeIsConsistent(fn) );

    if(DEBUG_MODE) {
      for(fn_ptr child = node->first_child; child != container.end() 
	    && child->parent == fn; ++child) {
	assert(nodeIsConsistent(child));
      }
    }

    // Now send all excess mass to the root, as much as possible.
    // Things are not made to reflect the node 

    fn->flow_to_parent += 0;
    
    
  }    

  void makeConsistent(fn_ptr fn) {
    assert(! nodeIsConsistent(fn) );
    
    // Ensure all the children are consistent; otherwise this is pointless.
    if(DEBUG_MODE) {
      if(flow_outward) {
	for(fn_ptr child = node->first_child; child != container.end() && child->parent == fn; ++child) {
	  assert(nodeIsConsistent(child));
	}
      } else {
	assert(nodeIsConsistent(fn->parent));
      }
    }
    
    // Now go and reconcile everything
    if(flow_outward) {
      

    }
  }

  // Each time we hit a threshhold quantity of excess, send it back
  // down to the root in order to see 


  template <int direction> 
  inline void pullInNeighbor(fn_ptr parent, int dim) {
    const node_ptr node = fnode->node;

    node_ptr nn = lattice.neighbor(node, dim, direction);

    if(nn->on() == partition) {

      ftype pull_capacity = flowCapacity<partition, direction>(nn, node, dim);

      if(pull_capacity > 0) {

	const bool in_needy_neighborhood = (node->need() > 0);
      
	if(nn->need() > 0) {
	  fn_ptr = fnc.newFlowNode(nn, parent, DirDim(dim, direction) );
	  fn_ptr.need = fnode->need() + need;

	} else if (nn->excess() > 0) {
	
	
	}
      }


      ftype required_need = fnode->required_need;
      ftype desired_need = fnode->desired_need + node->need();

      ftype need = required_need + desired_need;


    

  
  
      && !nn->isVisited() && nn->fowardPullCapacity(i) != 0) {
  
      // If we're in a need

  
    }


  };


  ////////////////////////////////////////////////////////////////////////////////
  // All the data structures for the active node reconcilliation


  // Returns true if this saturated that edge, which means we can
  // ignore this node now...

  bool pullAvailableAlongNonTreePath(fn_ptr parent, dn_ptr pull_dn) {

    node_ptr pull_node = pull_dn->node;
    node_ptr parent_node = parent->node;

    ftype amount = min(flowCapacity<partition>(pull_node, parent_node, pull_dn->step_to_parent),
			pull_node->visitor()->mass);

    pull_node->visitor()->mass -= amount;
    pullExcess(parent_node, pull_node, pull_dn->step_to_parent, amount);
    parent_node->mass += amount;
  }

  inline void pullAvailableAlongNonTreePath(fn_ptr parent, dn_ptr pull_dn) {

    node_ptr pull_node = pull_dn->node;
    node_ptr parent_node = parent->node;

    ftype amount = min(flowCapacity<partition>(pull_node, parent_node, pull_dn->step_to_parent),
			pull_node->visitor()->mass);

    pull_node->visitor()->mass -= amount;
    pullExcess(parent_node, pull_node, pull_dn->step_to_parent, amount);
    parent_node->mass += amount;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Movement operations
  
  template <int direction>
  inline bool _en_handleVisitedNode(fn_ptr fn, node_ptr nn, int dim) {
    // Currently, do nothing
    return false;
  }

  template <int direction>
  inline bool _en_preview_handleNewNode(fn_ptr fn, node_ptr nn, int dim) {
    // Attempt to pull in mass from these nodes; if it saturates the
    // edges, great.

    if(nn->need() > 0) {
      
      if(fn->node->need() == 0) {
	fn->external_need_at_node_and_above += 
	  min(nn->need(), e.flowCapacity<partition, direction>(fn->node, nn, dim));
	return false;
      } else {
	return true;
      }
    } else {
      ftype excess = nn->excess<partition>();
      ftype max_flow = e.flowCapacity<partition, -direction>(nn, fn->node, dim);

      bool is_saturated = (max_flow <= excess);

      ftype pull_amount = is_saturated ? max_flow : excess;
			     
      e.pullExcess<partition, direction>(fn->node, nn, dim, pull_amount);
      fn->mass += pull_amount;

      return !is_saturated; 
    }
  }

  template <int direction> 
  inline void _en_handleVerifiedNode(fn_ptr fn, node_ptr nn, int dim) {
    if(nn->need() > 0) {
      // Here we 
      
    }
  }

  void expandNode(fn_ptr fn) {
      
    // First deal with the nodes that we can actually pull things
    // from.  It may be that these nodes are 

    Array<int, 2*n_dimensions> node_reexamine(false);

    // First, we go and pull in what we can from neighbors.  It's nice
    // if we can saturate edges; then we don't have to worry about
    // those nodes.

    for(int dim = 0; dim < n_dimensions; ++dim) {
      node_ptr nn = lattice.neighbor(node, dim, -1);
      
      if(nn->on() == partition) { 
	if(nn->isVisited())
	  node_reexamine[dim] = _en_handleVisitedNode<-1>(fn, nn, dim);
	else
	  node_reexamine[dim] =_en_preview_handleNewNode<-1>(fn, nn, dim);
      }
    }

    for(int dim = 0; dim < n_dimensions; ++dim) {
      node_ptr nn = lattice.neighbor(node, dim, 1);
      
      if(nn->on() == partition) { 
	if(nn->isVisited())
	  node_reexamine[n_dimensions + dim] = _en_handleVisitedNode<1>(fn, nn, dim); 
	else
	  node_reexamine[n_dimensions +dim] =_en_preview_handleNewNode<1>(fn, nn, dim);
      }
    }

    // Now go through and 

	ftype capacity = 


	  ftype pull_amount = min(nn->forwardPullCapacity(i), nn->excess());
	  
	if(pull_amount > 0) {
	  pullExcess<partition, -1>(nn, node, min(pull_amount, desired_need) );
	} else if { 
	}
      }
    }
  }
    


 

  }

  ////////////////////////////////////////////////////////////////////////////////
  // High-level operations

  void pullEverythingToCenter() {
    
    _debugCheckContainerConsistency();

    // First pull out all the 

    unsigned int current_depth = container.back().depth;

    vector<list<fn_ptr> > queue(current_depth + 1);

    // Fill out 
    for(fn_ptr fn = tree_tips.begin(); fn != tree_tips.end(); ++fn) {
      queue[fn->depth].push_back(fn);
    }

    {
      fn_ptr fn = container.end(); --fn;
      
      for(; fn.depth == current_depth; --fn)
	queue[current_depth].push_back(fn);
    }
    
    for(; current_depth != 0; --current_depth) {
      
      for(fn_ptr fn = queue[current_depth].begin(); fn != queue[current_depth].end(); ++fn) {

	// Send as much mass as possible to the parent
	ftype inward_flow_cap = (flowCapacity<partition>(fn->node, fn->parent->node, 
							  fn->step_to_parent) 
				  - fn->flow_to_parent);

	ftype push_amount = min(fn->mass, forward_flow_cap);

	fn->mass -= push_amount;
	fn->flow_to_parent += push_amount;
	fn->parent->mass += push_amount;

	// This requires that the parent's need_above_node has been
	// set to zero by a previous process, either creation or the
	// out-propegation.
	ftype outward_flow_cap = (flowCapacity<partition>(fn->parent->node, fn->node, 
							   fn->step_to_parent)
				   + fn->flow_to_parent);
	
	fn->parent->need_at_node_and_above += min(fn->need_at_node_and_above, outward_flow_cap);
	fn->parent->need_at_branch_tips += min(fn->need_at_branch_tips, outward_flow_cap);
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Debug routines

  void _debugCheckContainerConsistency() {
    unsigned int cur_depth = 0;

    for(fn_ptr it = container.begin(); it != container.end(); ++it) {
      if(it->depth != cur_depth) {
	assert_equal(cur_depth+1, it->depth);
	++cur_depth;
      }
    }

    for(fn_ptr it = tree_tips.begin(); fn != tree_tips.end(); ++fn) {
      assert(it->first_child == container.end());
    }

  }

  void _debugCheckFNodeConsistency(fn_ptr fn) {
    // Ensure that the mass at that node is equal to the excess at
    // that node plus all the mass coming in, minus that leaving to
    // the parents

    _debugCheckNodeConsistency(fn->node);

    ftype mass = fn->node->excess();
    
    for(fn_ptr child = fn.first_child; child != container.end() && child.parent == fn; ++child) {
      mass += child->flow_to_parent;
    }
    
    mass -= fn->flow_to_parent;

    assert_equal(mass, fn->mass);
  }

////////////////////////////////////////////////////////////////////////////////
// CRAP code

  template <int partition>
  void reconcileActiveNode(node_ptr node) {
      
    if(partition_flag)
      assert( partition_flag ? (node->reduction < 0) : (node->reduction > 0) );

    // Set up the first one 
    FlowNodeContainer flow_nodes;
      
    flow_nodes.push_back(FlowNode(node));

    flow_nodes.back().required_need = abs(node->reduction);
    flow_nodes.back().capacity_to_base = need.required_need;

    node->visitor = &(flow_nodes.back());

    appendLocalActiveNodes();
  }


};

