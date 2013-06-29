#ifndef _TV_SOLVER_H_
#define _TV_SOLVER_H_

#include "energy.hpp"
#include "kernels/kernels.hpp"
#include <cmath>

namespace latticeQBP {

  static inline long toLong(double x) { return long(round(x*(1 << 12))); }
  static inline double toDbl(long x) { return double(x) *(1 / (1 << 12)); }

  template <typename Kernel, typename dtype = long>
  vector<double> calculate2dTV(size_t nx, size_t ny, 
                               double *function, double lambda) {

    static_assert(Kernel::is_geocut_applicable,
                  "Kernel is not valid for GeoCuts or TV Minimization.");
    static_assert(Kernel::n_dimensions == 2, "Currently only dealing with 2d stuff.");

    typedef LatticeLevelReductions<2, Kernel, dtype> rsolver_type;
    typedef typename rsolver_type::index_vect index_vect;

    rsolver_type rsolver(index_vect({nx, ny}));

    for(auto pwfi = rsolver.getPairwiseFillingIterator(); !pwfi.done(); ++pwfi) {
        
      size_t src_idx_x = pwfi.latticeCoordOf1()[0];
      size_t src_idx_y = pwfi.latticeCoordOf1()[1];
      size_t idx_src = nx*src_idx_y + src_idx_x;

      size_t dest_idx_x = pwfi.latticeCoordOf2()[0];
      size_t dest_idx_y = pwfi.latticeCoordOf2()[1];
      size_t idx_dest = nx*dest_idx_y + dest_idx_x;

      long pwf = toLong(0.5*lambda*pwfi.geocutEdgeWeight() 
                        * abs(function[idx_src] - function[idx_dest]));

      pwfi.addPairwisePotential(0, pwf, pwf, 0);
      pwfi.addUnaryPotentialTo1(0, function[idx_dest]);
    }

    rsolver.run();

    vector<double> res(nx * ny);

    size_t i = 0;
    for(auto it = rsolver.getLattice().indexIterator(); !it.done(); ++it) {
      res[i++] = toDbl(rsolver.getLattice()(it.coords())->reduction);
    }

    return res;
  }
};

#endif /* _TV_SOLVER_H_ */
