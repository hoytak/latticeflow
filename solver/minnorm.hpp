
#ifndef _ENERGY_PRIVATE_MEMBER_
#error "File minnorm.hpp cannot be included directly."
#else

////////////////////////////////////////////////////////////////////////////////
// Run 
 
void runPRMinNorm() {

  greedy();

  cout << "Size of node = " << sizeof(Node) << " bytes." << endl;

  // Set up the flow algorithms
  PRFlow<0> pr_flow_0(*this);
  PRFlow<1> pr_flow_1(*this);

  // Set up the processing queue
  // deque<node_ptr> node_queue;

  for(size_t iteration = 0; iteration < 1000; ++iteration) {
    
    cout << "Iteration " << (iteration + 1) << ": " << endl;

    bool r1 = pr_flow_1.runFullPass();
    bool r0 = pr_flow_0.runFullPass();

    if(r1 && r0) {
      return;
    }
  }
}

void runPRHeapMinNorm() {

  greedy();

  // Set up the flow algorithms
  PRFlow<0> pr_flow_0(*this);
  PRFlow<1> pr_flow_1(*this);

  // Set up the processing queue
  // deque<node_ptr> node_queue;

  size_t flip_count_0 = 0, flip_count_1 = 0;

  for(size_t iteration = 0; iteration < 100000; ++iteration) {
    
    dtype base_reduction = 0;
    node_ptr max_r_n = lattice.end();

    // Find the one that is most violating things 
    for(node_ptr n = lattice.begin(); n != lattice.end(); ++n) {
      dtype r = (n->on() ? 1 : -1) * n->reduction;
      
      if(r > base_reduction) {
	base_reduction = r;
	max_r_n = n;
      }
    }

    if(max_r_n == lattice.end()) 
      break;

    if(max_r_n->on())
      flip_count_1 += pr_flow_1.run(max_r_n);
    else 
      flip_count_0 += pr_flow_0.run(max_r_n);
  }

  cout << "Running, did " << flip_count_0 << " 0 -> 1 flips and " << flip_count_1 << " 1 -> 0 flips." << endl;

}

void runTabuMinNorm() {

  TimeTracker tt;

  tt.start();

  TabuSearch ts(*this);

  size_t active_nodes = 0;
  
  for(node_ptr n = lattice.begin(); n != lattice.end(); ++n) {
    if(n->on() ? (n->reduction > 0) : (n->reduction < 0) ) {
      ++active_nodes;
      ts.addActiveNode(n);
    }
  }

  tt.stop();
  cout << "Took " << tt.asString() << " to set up tabu search with " 
       << active_nodes << " active nodes." << endl;

  tt.reset();
  tt.start();

  ts.run();

  tt.stop();
  cout << "Took " << tt.asString() << " to run tabu search. " << endl;

  tt.reset();
  tt.start();

  // Set up the flow algorithms
  PRFlow<0> pr_flow_0(*this);
  PRFlow<1> pr_flow_1(*this);

  // Set up the processing queue
  // deque<node_ptr> node_queue;

  for(size_t iteration = 0; iteration < 1000; ++iteration) {
    
    cout << "Iteration " << (iteration + 1) << ": " << endl;

    bool r1 = pr_flow_1.runFullPass();
    bool r0 = pr_flow_0.runFullPass();

    if(r1 && r0) {
      return;
    }
  }

  tt.stop();
  cout << "Took " << tt.asString() << " to run PR flow. " << endl;

}
#endif
