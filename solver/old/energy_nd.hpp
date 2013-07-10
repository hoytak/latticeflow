#ifndef _ENERGY2D_H_
#define _ENERGY2D_H_

#include "energy.hpp"

using latticeQBP::LatticeEnergy;

template <typename dtype, typename ftype> 
class Energy1D : public LatticeEnergy<1, dtype, ftype> {
public:
  Energy1D(size_t n_x) 
    : LatticeEnergy<1,dtype,ftype>(make_array_1(n_x) )
  {}

  void addE1(size_t xi, dtype E0, dtype E1) {
    LatticeEnergy<1,dtype,ftype>::addE1(make_array_1(xi), E0, E1);
  }
  
  // Indexed by the lower node
  void addE2(size_t xi, dtype E00, dtype E01, dtype E10, dtype E11) {
    LatticeEnergy<1,dtype,ftype>::addE2(make_array_1(xi), 0, E00, E01, E10, E11);
  }

  bool on(size_t xi) const {
    return LatticeEnergy<1,dtype,ftype>::on(make_array_1(xi));
  }
};

template <typename dtype, typename ftype> 
class Energy2D : public LatticeEnergy<2, dtype, ftype> {
public:
  Energy2D(size_t n_x, size_t n_y) 
    : LatticeEnergy<2,dtype,ftype>(make_array_2(n_x, n_y) )
  {}

  void addE1(size_t n_x, size_t n_y, dtype E0, dtype E1) {
    addE1(make_array_2(n_x, n_y), E0, E1);
  }

  void addE2(size_t n_x, size_t n_y, dtype E00, dtype E01, dtype E10, dtype E11) {
    addE2(make_array_2(n_x, n_y), E00, E01, E10, E11);
  }
};

template <typename dtype, typename ftype> 
class Energy3D : public latticeQBP::LatticeEnergy<3, dtype, ftype> {
public:
  Energy3D(size_t n_x, size_t n_y, size_t n_z) 
    : LatticeEnergy<3,dtype,ftype>(make_array_3(n_x, n_y, n_z) )
  {}

  void addE1(size_t n_x, size_t n_y, size_t n_z, dtype E0, dtype E1) {
    addE1(make_array_3(n_x, n_y, n_z), E0, E1);
  }

  void addE2(size_t n_x, size_t n_y, size_t n_z, 
	     dtype E00, dtype E01, dtype E10, dtype E11) {
    addE2(make_array_3(n_x, n_y, n_z), E00, E01, E10, E11);
  }
};

#endif /* _ENERGY2D_H_ */
