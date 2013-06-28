
#ifndef _ENERGY_PRIVATE_MEMBER_
#error "File nodes.hpp cannot be included directly."
#else

// Very simple
void greedy() {

  TimeTracker tt;

  tt.start();

  // Could be optimized a LOT
  bool changed;
  do{
    changed = false;
    for(node_ptr n = lattice.begin(); n != lattice.end(); ++n) {
      if(n->gain > 0) {
	flipNode<3>(n);
	changed = true;
      }
    }
  }while(changed);

  tt.stop();

  cout << "Greedy finished in " << tt.asString() << "." << endl;
}


#endif
