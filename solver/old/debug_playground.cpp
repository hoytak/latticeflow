#ifdef ENABLE_OPENMP
#undef ENABLE_OPENMP
#endif

#include "energy_nd.hpp"

int main(int argn, char** argv)
{ 
  Energy1D<int, double> e(3);
  
  e.addE1(0, 0, 10);
  e.addE1(1, 0, 10); 
  e.addE1(2, 0, 35);
  e.addE2(0, 0,0,0,-30);
  e.addE2(1, 0,0,0,-30);

  e.run();

  assert(e.on(0));
  assert(e.on(1));
  assert(!e.on(2));
}


