#ifndef _TV_FUNCTION_MAP_H_
#define _TV_FUNCTION_MAP_H_

#include "../lattices/lattice.hpp"

namespace latticeQBP {
  
  using namespace std;  
  
  template <typename dtype, int n_dimensions>
  class FunctionMap : LatticeArray<dtype, n_dimensions> {

    typedef LatticeArray<dtype, n_dimensions> Base;
    typedef Base::index_vect index_vect;

    FunctionMap(index_vect dimensions) 

  };
};

#endif /* _TV_FUNCTION_MAP_H_ */
