// File meant to be included only from within the energy struct

#ifndef _ENERGY_PRIVATE_MEMBER_
#error "File energy_debug.hpp cannot be included directly."
#else

void _debugPrintState() const {
  cout << "State (value=" << current_value << "): ";
      
  for(node_ptr n = lattice.begin(); n != lattice.end(); ++n) 
    cout << (n->on() ? '+' : '-');
      
  cout << endl;
}

////////////////////////////////////////////////////////////////////////////////
// Functions for doing the bookkeeping with the nodes, setting up
// the flips, etc.
void _debugVerifyNodeConsistency(node_ptr node) { 

#ifndef NDEBUG

  if(!(node >= lattice.begin() && node < lattice.end()))
    return;

  assert_lt(node->gain, Node::pinNodeConstant() / 2);

  bool pinned = (node->isPinned());
  
  if(pinned)
    node->unpin();

  // Make sure the gain is what we'd like
      
  dtype value_on = node->q;

  for(size_t i = 0; i < n_dimensions; ++i) {
    {
      node_ptr nn = lattice.neighbor(node, i, -1);
      if(lattice.isValidNode(nn) && nn->on())
	value_on += nn->edges[i];
    }

    {
      node_ptr nn = lattice.neighbor(node, i, 1);
      if(lattice.isValidNode(nn) && nn->on())
	value_on += node->edges[i];
    }
  }

  if(node->on()) {
    assert_equal(node->gain, value_on);
  } else {
    assert_equal(node->gain, -value_on);
  }

  if(pinned)
    node->pin();

#endif
}

void _debugVerifyNodeReduction(node_ptr node, bool print_always = false) { 
#ifndef NDEBUG

  if(!(node >= lattice.begin() && node < lattice.end()))
    return;

  ftype r_value = 2*node->q;

  if(node->on()) {
	
    for(size_t i = 0; i < n_dimensions; ++i) {
      {	  
	node_ptr nn = lattice.neighbor(node, i, -1);
	    
	assert_leq(abs(nn->edges[i]), 100000);
	assert_leq(abs(nn->alpha[i]), 100000);
	    
	assert_leq(abs(nn->alpha[i]), abs(nn->edges[i]));

	if(nn->on()) {
	  r_value += nn->edges[i];
	  r_value += nn->alpha[i];
	}
      }

      {
	node_ptr nn = lattice.neighbor(node, i, 1);

	assert_leq(abs(node->edges[i]), 100000);
	assert_leq(abs(node->alpha[i]), 100000);

	assert_leq(abs(node->alpha[i]), abs(node->edges[i]));

	if(nn->on()) {
	  r_value += node->edges[i];
	  r_value -= node->alpha[i];
	}
      }
    }

  } else {

    for(size_t i = 0; i < n_dimensions; ++i) {
      {	  
	node_ptr nn = lattice.neighbor(node, i, -1);

	assert_leq(abs(nn->edges[i]), 100000);
	assert_leq(abs(nn->alpha[i]), 100000);

	assert_leq(abs(nn->alpha[i]), abs(nn->edges[i]));

	if(nn->on()){
	  r_value += 2*nn->edges[i];
	} else {
	  r_value += nn->edges[i];
	  r_value += nn->alpha[i];
	}
      }

      {
	node_ptr nn = lattice.neighbor(node, i, 1);

	assert_leq(abs(node->edges[i]), 100000);
	assert_leq(abs(node->alpha[i]), 100000);

	assert_leq(abs(node->alpha[i]), abs(node->edges[i]));

	if(nn->on()){
	  r_value += 2*node->edges[i];
	} else {
	  r_value += node->edges[i];
	  r_value -= node->alpha[i];
	}
      }
    }
  }

  bool _node_error = abs(node->reduction - r_value) > 1e-4;

  if(_node_error)
    cout << ">>>> REDUCTION BOOKKEEPING ERROR <<<< " << endl;


  if(_node_error || print_always) {
	
    cout << "Node index = " << (node - lattice.begin()) << endl;
    cout << "r_value = " << r_value << "; node->reduction = " << node->reduction << endl;

    cout << "(" << (node - lattice.begin())
	 << "); state = " << node->state
	 << "; q = " << node->q
	 << ": reduction = " << node->reduction 
	 << endl;

    for(size_t i = 0; i < n_dimensions; ++i) {
	  
      node_ptr nn = lattice.neighbor(node, i, -1);

      cout << "  -" << i 
	   << "; state = " << nn->state
	   << "; q = " << nn->q
	   << "; edge = " << nn->edges[i]
	   << "; alpha = " << nn->alpha[i]
	   << ": reduction = " << nn->reduction 
	   << endl;
    }

    for(size_t i = 0; i < n_dimensions; ++i) {
	  
      node_ptr nn = lattice.neighbor(node, i, 1);

      cout << "  +" << i 
	   << "; state = " << nn->state
	   << "; q = " << nn->q
	   << "; edge = " << node->edges[i]
	   << "; alpha = " << node->alpha[i]
	   << ": reduction = " << nn->reduction 
	   << endl;
    } 
	
    if(_node_error)
      abort();
  }
#endif
}
#endif
