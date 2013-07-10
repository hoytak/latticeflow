#ifndef _ENERGY_PRIVATE_MEMBER_
#error "File tabusearch.hpp cannot be included directly."
#else

#include "greedy_heap.hpp"

#define TABU_PENALTY_EXPFACTOR (size_t(5))
#define TABU_PENALTY_MINIMUM   (size_t(16))

class TabuSearch {
public:
  
  TabuSearch(LatticeEnergy& _e)
    : e(_e)
    , lattice(e.lattice)
    , n_active(0)
    , state_tracker(e.lattice.size(), e.lattice.size() / 10)
    , p_gain(0)
    , fb(0)
  {
  }
  
private:
  
  LatticeEnergy& e;
  Lattice& lattice;

  size_t n_active;
  GreedyHeap run_heap;

  deque<node_ptr> tabu_ring;
  typedef typename deque<node_ptr>::iterator tabu_ring_iterator;

  StateTracker<unsigned int> state_tracker;
  
  dtype p_gain;
  long prev_best_value;

  FastBernoulli fb;

  ////////////////////////////////////////////////////////////////////////////////

  node_ptr top() {
    node_ptr n;
    run_heap.top(n, true);
    return n;
  }
  
  void pop() {
    run_heap.pop();
  }

  size_t tabuPenalty() const {
    size_t v = max(min( (n_active >> 1), TABU_PENALTY_MINIMUM), (n_active >> TABU_PENALTY_EXPFACTOR));
    return v;
  }

  void randomJump(size_t np) {
    clearTabuRing();

    vector<node_ptr> perturb_list(np);

    for(size_t i = 0; i < np; ++i) {
      bool done = run_heap.top(perturb_list[i]);

      if(unlikely(done)) {
	perturb_list.resize(i);
	break;
      }
    }

    for(size_t i = 0; i < perturb_list.size(); ++i) {
      if(fb())
	flipNode(perturb_list[i], false);
      else
	run_heap.push(perturb_list[i]);
    }
    
    runGreedy();
  }

  inline void flipNode(node_ptr n, bool pin) {

    assert(!n->isPinned());

    p_gain = n->gain;

    e.template flipNode<3>(n);

    assert(!n->isPinned());

    state_tracker.flip(n - lattice.begin());
    
    for(int i = 0; i < 2*n_dimensions; ++i) {
      node_ptr nn = lattice.neighbor(n, i);
      if(e.isConnected(n, nn, i)) {
	assert(nn >= lattice.begin() && nn < lattice.end());

	if(! nn->visited) {
	  nn->visited = true;
	  ++n_active;
	}
	
	if(!nn->isPinned())
	  run_heap.push(nn);
      }
    }
    
    if(pin) {
      n->pin();
      tabu_ring.push_back(n);
    } else {
      run_heap.push(n);
    }
  }
  
  void unpinNode(node_ptr n) {
    // cout << "unpinned node " << (n - lattice.begin()) << "; " << "freeze_ring size = " << tabu_ring.size() << endl;
    n->unpin();
    run_heap.push(n);
  }
  
  void advanceFreezeRing() {
    for(int i = 0; i < 3; ++i) {
      if(tabu_ring.size() > tabuPenalty()) {

 	node_ptr n = tabu_ring.front();
	tabu_ring.pop_front();
	unpinNode(n);
      }
    }
  }

  void clearTabuRing() {
    for(tabu_ring_iterator it = tabu_ring.begin(); it != tabu_ring.end(); ++it)
      unpinNode(*it);

    if(DEBUG_MODE) {
      for(node_ptr n = lattice.begin(); n != lattice.end(); ++n)
	assert(!n->isPinned());
    }

    tabu_ring.clear();
  }

  void runGreedy() {
    clearTabuRing();
    
    node_ptr n;

    while(true) {
      run_heap.top(n, true);
      if(likely(n->gain > 0)) {
	run_heap.pop();
	flipNode(n, false);
      } else {
	break;
      }
    }
  }

  bool atBend(node_ptr next_node) const {
    return next_node->gain < 0 && p_gain >= 0; 
  }

  bool hitNewBestValue() const {
    return e.current_value < prev_best_value;
  }

  void descendToNewLocalOptimal() {
    assert_lt(e.current_value, prev_best_value);
    runGreedy();
    assert_lt(e.current_value, prev_best_value);
    
    prev_best_value = e.current_value;
  }

  bool atKnownBest() const {
    return e.current_value == prev_best_value;
  }

public:

  void addActiveNode(node_ptr n) {

    e._debugVerifyNodeConsistency(n);

    run_heap.push(n);
    ++n_active;
  }
  
  void run() {

    size_t iteration = 0;

    while(true) {

    next_iteration:

      advanceFreezeRing();

      ++iteration;

      if(unlikely(iteration == 10*lattice.size())) {
	runGreedy();
	goto done;
      }

      // Get the top node, then see what to do with it.
      node_ptr n = top();

      if(atBend(n)) {
	
	// Only place the state tracker should register something.

	unsigned int& value = state_tracker.value();
	
	switch(++value) {
	case 1:
	  break;
	case 2:
	  cout << (atKnownBest() ? "<B2>" : "2");
	  pop();
	  flipNode(n, false);
	  randomJump((atKnownBest() ? 4 : 2) * tabuPenalty());
	  goto next_iteration;
	case 3:
	  cout << (atKnownBest() ? "<B3 Done>" : "3");

	  if(unlikely(atKnownBest()))
	    goto done;

	  flipNode(n, false);
	  randomJump(4*tabuPenalty());
	  goto next_iteration;

	default:
	  cout << "4";

	  assert(!atKnownBest());

	  randomJump(n_active / 2);
	  goto next_iteration;
	}
      }

      // Okay, well, this is good!
      flipNode(n, true);

      // If the current value is better than any previously known,
      // then run greedy, record the solution, and then move on from
      // there.
      if(hitNewBestValue()) {
	descendToNewLocalOptimal();
	goto next_iteration;
      }
    }

  done:
    
    cout << endl;
    
    for(node_ptr n = lattice.begin(); n != lattice.end(); ++n) {
      n->visited = false;
    }
  }

};

#endif
