#ifndef _TV_SOLVER_H_
#define _TV_SOLVER_H_

#include "../kernels/kernels.hpp"
#include "../common.hpp"
#include "../lattices/lattice.hpp"
#include "../parametricflow/mn_netflow_solver.hpp"
#include "tv_regpath.hpp"
#include "tv_push_relabel.hpp"

#include <cmath>
#include <iostream>
#include <cstdint>
#include <memory>

namespace latticeQBP {
  
  using namespace std;

  template <typename Kernel, typename dtype = long> 
  class TVSolver {
  public:
    static_assert(Kernel::is_geocut_applicable,
                  "Kernel is not valid for GeoCuts or TV Minimization.");
    static_assert(Kernel::n_dimensions == 2, "Currently only dealing with 2d stuff.");

    static_assert(sizeof(dtype) >= 32, "This class does not work with lower precision data type; use dtype of 32 bits or more.");

    
    // Set up the kernel stuff. 
    static constexpr int nf_parameters = (NF_ON_BY_REDUCTION 
                                          | NF_ENABLE_KEY_PARTIONING 
                                          | NF_ADJUSTMENT_MODE);


    typedef TVFlowNode<Kernel, dtype> Node;
    typedef KernelLattice<Node, Kernel::n_dimensions, Kernel> Lattice;
    typedef typename Lattice::index_vect index_vect;
    typedef typename Lattice::value_ptr node_ptr; 
    typedef typename CompType<dtype>::Type comp_type;
    typedef TV_PRFlow<dtype, Lattice> PRSolver;
    typedef TVRegPathManager<dtype, PRSolver> TVRegPathManager;

    
  private:

    ////////////////////////////////////////////////////////////////////////////////
    // General information
    size_t nx, ny;
    double *function;
    Lattice lattice;
    dtype max_lambda;
    TVRegPathManager reg_path;

  public:

    TVSolver(size_t _nx, size_t _ny, double *_function) 
      : nx(_nx)
      , ny(_ny)
      , function(_function)
      , lattice(index_vect({nx, ny}))
      , reg_path(lattice)
    {
      ////////////////////////////////////////////////////////////
      // Initialize the lattice 

      for(auto it = lattice.indexIterator(); !it.done(); ++it) {
        const index_vect& idx = it.coords();
        lattice(idx)->setBaseFunctionValue(function[nx*idx[0] + idx[1]]);
      }
      
      Node::initTVLattice(lattice);
    }


  private:
    ////////////////////////////////////////////////////////////////////////////////
    // The runners 


    PRSolver pr_solver;

    ////////////////////////////////////////////////////////////////////////////////
    // Initializing all the sets 

    struct NodePivot;
    typedef shared_ptr<set<node_ptr> > nodeset_ptr;

    struct RegPathSegment {

    ////////////////////////////////////////////////////////////
    // This is simply a storage container for node pivots.  
    list<RegPathSegment> node_pivot_stack;
    
    RegPathSegment* getNewRegPathSegment() {
      node_pivot_stack.push_back(RegPathSegment());
      return &node_pivot_stack.back();
    }

    ////////////////////////////////////////////////////////////
    // The map for the path of this node starting at max_lambda
    vector<RegPathSegment*> start_node_map;

    ////////////////////////////////////////////////////////////
    // The queue for building these things.  Each of these is 

    struct PivotPoint {
      dtype event_lambda; 
      bool operator<(const PivotPoint& p) const {
        return event_lambda < p.event_lambda;
      }
    };
    
    struct SplitPivotPoint : public PivotPoint {
      RegPathSegment* pivot_node_spliting;
      shared_ptr<vector<node_ptr> > split_1_ptr, split_2_ptr;
    };
    
    struct JoinPivotPoint : public PivotPoint { 
      RegPathSegment *rps_1, *rps_2;
    };

    priority_queue<SplitPivotPoint> splits;
    priority_queue<JoinPivotPoint> joins;

    list<RegPathSegment*> active_nodes;
    
    ////////////////////////////////////////
    // Functions to test things. 
    
                     
    ////////////////////////////////////////////////////////////////////////////////
    // Information from the 

  };

  ////////////////////////////////////////////////////////////////////////////////
  // Get a convenience function for the simple 2d case

  template <typename Kernel, typename dtype = long>
  vector<double> calculate2dTV(size_t nx, size_t ny, 
                               double *function, double lambda) {

    static_assert(Kernel::is_geocut_applicable,
                  "Kernel is not valid for GeoCuts or TV Minimization.");
    static_assert(Kernel::n_dimensions == 2, "Currently only dealing with 2d stuff.");

    typedef LatticeLevelReductions<2, Kernel, dtype> rsolver_type;
    typedef typename rsolver_type::index_vect index_vect;

    // cout << "nx = " << nx << "; ny = " << ny << endl;

    double min_x = *min_element(function, function + nx*ny);
    double max_x = *max_element(function, function + nx*ny);

    double w = max_x - min_x;

    const double conversion_factor = 
      (double(dtype(1) << (sizeof(dtype)*8 - 16))
       / (max(1.0, lambda) * (max_x - min_x)));

    auto toDtype = [conversion_factor](double x) {
      return dtype(round(x * conversion_factor));
    }; 

    auto toDbl = [conversion_factor](dtype x) {
      return double(x) / conversion_factor;
    }; 

    rsolver_type rsolver(index_vect({ny, nx}));

    for(auto ufi = rsolver.getUnaryFillingIterator(); !ufi.done(); ++ufi) {

      size_t idx_y = ufi.latticeCoord()[0];
      size_t idx_x = ufi.latticeCoord()[1];
      size_t idx = nx*idx_y + idx_x;

      dtype fv = toDtype(lambda*function[idx]);

      ufi.addUnaryPotential(0, fv);
    }

    for(auto pwfi = rsolver.getPairwiseFillingIterator(); !pwfi.done(); ++pwfi) {
        
      size_t src_idx_y = pwfi.latticeCoordOf1()[0];
      size_t src_idx_x = pwfi.latticeCoordOf1()[1];
      size_t idx_src = nx*src_idx_y + src_idx_x;

      size_t dest_idx_y = pwfi.latticeCoordOf2()[0];
      size_t dest_idx_x = pwfi.latticeCoordOf2()[1];
      size_t idx_dest = nx*dest_idx_y + dest_idx_x;

      dtype pwf = toDtype(0.5*pwfi.geocutEdgeWeight() 
                        * abs(function[idx_src] - function[idx_dest]));

      pwfi.addPairwisePotential(0, pwf, pwf, 0);

      // cout << pwfi.latticeCoordOf1() << " -> " << pwfi.latticeCoordOf2() 
      //      << ": pwf = " << pwf << endl;
    }

    rsolver.run();

    vector<double> res(nx * ny);

    size_t i = 0;
    for(auto it = rsolver.getLattice().indexIterator(); !it.done(); ++it) {
      // cout << it.coords() << ": r = " << rsolver.getLattice()(it.coords())->level() << endl;
      res[i++] = toDbl(rsolver.getLattice()(it.coords())->level()) 
        / (1e-32 + lambda);
    }

    return res;
  }
};

#endif /* _TV_SOLVER_H_ */
