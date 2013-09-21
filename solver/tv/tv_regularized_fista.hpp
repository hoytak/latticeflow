#ifndef _TV_REGULARIZED_FISTA_H_
#define _TV_REGULARIZED_FISTA_H_

#define _DEBUG_IN_TV_SOLVER_INCLUDES_

#include "../common.hpp"
#include "tv_push_relabel.hpp"
#include "tv_regpath_node.hpp"
#include "../lattices/lattice.hpp"
#include "../parametricflow/pf_solver.hpp"
#include "../parametricflow/pf_policies.hpp"

#undef _DEBUG_IN_TV_SOLVER_INCLUDES_

#include <set>
#include <list>
#include <vector>
#include <deque>
#include <queue>
#include <unordered_map>
#include <stdexcept>
#include <memory>

namespace latticeQBP { 
  namespace LFInternal {

    template <
      class Kernel,
      typename dtype,
      typename DataLattice,
      typename FunctionAndDerivative>
    class _TVRegularizedFista {
    public:

      typedef LatticeLevelReductions<Kernel::n_dimensions, Kernel, dtype> rsolver_type;
      typedef typename rsolver_type::index_vect index_vect;
      
      static_assert(Kernel::is_geocut_applicable,
                    "Kernel is not valid for GeoCuts or TV Minimization.");

      _TVRegularizedFista(DataLattice& value,
                          FunctionAndDerivative& _f) 
        : current_x(value)
        , f(_f)
      { }

      void run(double tv_reg_p, 
               double conv_threshhold,
               size_t max_iterations,
               double starting_L
               ) {
        
      }

    private:

      ////////////////////////////////////////////////////////////////////////////////
      // Data level functions to help with the solving


      // Does conversion stuff 
      double mid, conversion_factor;

      inline dtype data2dtype(double x) {
        return dtype(round(x * conversion_factor));
      }

      inline double dtype2data(dtype x) {
        return double(x) / conversion_factor;          
      }

      void setupConversionRoutines(const DataLattice& value, double tv_reg_p) {
        double min_x = *min_element(value.begin(), value.end());
        double max_x = *max_element(value.begin(), value.end());

        double w = max_x - min_x;
        mid = 0.5*(max_x + min_x);

        conversion_factor =
          (pow(2.0, (sizeof(dtype)*5))
           / (max(1.0, tv_reg_p) * (max_x - min_x)));
      }

      ////////////////////////////////////////////////////////////////////////////////
      // Stuff to move data back and forth. 

      void copyValueToSolver(const DataLattice& value) {
        
        


      }


      rsolver_type rsolver(value.shape());

    for(auto ufi = rsolver.getUnaryFillingIterator(); !ufi.done(); ++ufi) {

      size_t idx_x = ufi.coords()[0];
      size_t idx_y = ufi.coords()[1];

      dtype fv = toDtype(value[index_vect{idx_x, idx_y}] - mid);

      ufi.addUnaryPotential(0, fv / 2);
      assert_equal(2*fv, ufi.node()->r());
    }

    // for(auto it = rsolver.getLattice().vertexIterator(); !it.done(); ++it) {
    //   assert_almost_equal(function[it.nodeIndex()], mid + 0.5*toDbl(it->r()));
    // }

    for(auto pwfi = rsolver.getPairwiseFillingIterator(); !pwfi.done(); ++pwfi) {

      size_t src_idx_x = pwfi.coordsOf1()[0];
      size_t src_idx_y = pwfi.coordsOf1()[1];

      size_t dest_idx_x = pwfi.coordsOf2()[0];
      size_t dest_idx_y = pwfi.coordsOf2()[1];

      dtype pwf = toDtype(0.25*reg_p);

      pwfi.addPairwisePotential(0, pwf, pwf, 0);
    }

    rsolver.run();

    FuncMapPtr Rptr = FuncMapPtr(new LatticeArray<double, 2>(index_vect({nx, ny})));

    for(auto it = rsolver.getLattice().vertexIterator(); !it.done(); ++it) {
      (*Rptr)[it.coords()] = mid + 0.5*toDbl(it->r());
      // assert_almost_equal(function[it.nodeIndex()], (*Rptr)[it.coords()]);
    }

    return Rptr;
      DataLattice2d& current_x;
      
      FunctionAndDerivative& f;
    };
  }

  
  void runTVRegularizedFista(DataLattice2d& value,
                             FunctionAndDerivative& f,
                             double tv_reg_p, 
                             double conv_threshhold,
                             size_t max_iterations,
                             double starting_L
                             )
  {



  }
  }
}

#include "../common/debug_flymake_test.hpp"

#endif 
