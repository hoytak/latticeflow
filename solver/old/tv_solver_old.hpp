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

    
    
    
    typedef typename Lattice::index_vect index_vect;
    typedef typename Lattice::value_ptr node_ptr; 
    typedef typename CompType<dtype>::Type comp_type;
    
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
};

#endif /* _TV_SOLVER_H_ */
