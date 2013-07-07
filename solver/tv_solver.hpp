#ifndef _TV_SOLVER_H_
#define _TV_SOLVER_H_

#include "energy.hpp"
#include "kernels/kernels.hpp"
#include "tv_flow_node.hpp"
#include "numerical.hpp"

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

  private:

    ////////////////////////////////////////////////////////////////////////////////
    // General information
    size_t nx, ny;
    double *function;

    Lattice lattice;

    dtype max_lambda;

  public:

    TVSolver(size_t _nx, size_t _ny, double *_function) 
      : nx(_nx)
      , ny(_ny)
      , function(_function)
      , lattice(index_vect({nx, ny}))
      , pr_solver(lattice)
    {
      ////////////////////////////////////////////////////////////
      // Initialize the lattice 

      for(auto it = lattice.indexIterator(); !it.done(); ++it) {
        const index_vect& idx = it.coords();
        lattice(idx)->setBaseFunctionValue(function[nx*idx[0] + idx[1]]);
      }
      
      Node::initTVLattice(lattice);

      // Now we're ready to go...
    }


  private:
    ////////////////////////////////////////////////////////////////////////////////
    // The runners 

    typedef PRFlow<dtype, Lattice, 0, Node::nf_parameters> PRSolver;

    PRSolver pr_solver;

    ////////////////////////////////////////////////////////////////////////////////
    // Initializing all the sets 

    struct NodePivot;
    typedef shared_ptr<set<node_ptr> > nodeset_ptr;

    struct RegPathSegment {
      enum Mode {Terminal, Join, Split} rhs_mode, lhs_mode; 

      dtype rhs_lambda, lhs_lambda;
      dtype rhs_r, lhs_r;

      ////////////////////////////////////////////////////////////
      // If rhs_mode == Terminal:
      //    
      //    - nodeset gives the set of nodes along this path.  
      //    
      //    - rhs_nodes are both Null
      //    
      // 
      // If rhs_mode == Join:
      //    
      //    - nodeset is empty.  It is assumed to be a join of the two 
      //      nodes in rhs_nodes. 
      //    
      //
      // If rhs_mode == Split:
      //    
      //    - nodeset gives the nodes split off along this path.
      // 
      //    - If nodeset is empty, then rhs_lambda == lhs_lambda, and
      //      it's an intermediate structure of a 3+ way split.
      //    
      //    - rhs_nodes[0] gives the node that was split; rhs_nodes[1] is null.
      //    
      // 
      // Same with lhs, except two nodes for a split and one for a join. 
      //
      // Note that it is possible for lhs_lambda == rhs_lambda.  If
      // this happens, then the split / join is of more than one set
      // at the same lambda.  This is a rare occurance, so don't
      // optimize for it.
      //
      // While things are building, the nodeset contains the current
      // nodes for both the join and split cases.  But it's then
      // cleared out when the lhs is decided. 
      
      set<node_ptr> nodeset;
      
      // This is set up 
      Array<RegPathSegment*, 2> rhs_nodes;
      Array<RegPathSegment*, 2> lhs_nodes;
    };

    // This is simply a storage container for node pivots.  
    list<RegPathSegment> node_pivot_stack;

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
      RegPathSegment* split_point;
      shared_ptr<set<node_ptr> > split_1_ptr, split_2_ptr;
    };
      
    struct JoinPivotPoint : public PivotPoint { 
      RegPathSegment *rps_1, *rps_2;
    };

    priority_queue<SplitPivotPoint> splits;
    priority_queue<JoinPivotPoint> joins;

    list<RegPathSegment*> active_nodes;

    ////////////////////////////////////////////////////////////
    // Adding nodes in there 

    template <typename ForwardIterator> 
    void processNode(ForwardIterator it_start, ForwardIterator it_end, 
                     dtype lambda, dtype known_join_lambda = 0) {
      
      ////////////////////////////////////////////////////////////////////////////////
      // First, see if it's possible to solve the transhipment problem
      // at the lowest possible lambda

      StableAverage<dtype> fv_avg;
      
      for(ForwardIterator it = it_start; it != it_end; ++it) {
        node_ptr n = (*it);
        n->setFunctionValue(lattice, 0, known_join_lambda); 
        fv_avg.add(n->fv_predict());
      }
      
      dtype fv_offset = fv_avg.valueRoundedUp();

      for(ForwardIterator it = it_start; it != it_end; ++it) {      
        node_ptr n = (*it);        
        n->setOffset(lattice, fv_offset);
      }

      const uint key = 1021;

      pr_solver.prepareSection(it_start, it_end, key);

      pr_solver.runSection(it_start, it_end, key);

      // Now, see if there are any nodes that are on; if not, we are doing well. 
      

    }
      
                     
    template <typename ForwardIterator> 
    pair<dtype, dtype> getDTInformation(ForwardIterator it_start, ForwardIterator it_end) {
      
      // This is rather tricky, as we need to calculate a bunch of
      // averages without having numerical overflow issues. 

      StableAverage reductions, function_values;

      

    }

    ////////////////////////////////////////////////////////////////////////////////
    // Information from the 


    template <typename ForwardIterator> 
     get(ForwardIterator it_start, ForwardIterator it_end) {



    template <typename ForwardIterator> 
    void addStartingNode(ForwardIterator it_start, ForwardIterator it_end) {
      
      

    }





    ////////////////////////////////////////
    // Now for setting up the initial set.

    void fillIn() {


    }


    RegPathSegment* newInitialRegPathSegment(dtype lambda, 



    template <typename ForwardIterator> 
    void _addInitialSet


    void initSetsAtMaxLambda() { 
      
      deque<deque<node_ptr> > queue; 
      vector< 




    }

    // This calculates the TV for the 
    void run(double _max_lambda) {

      max_lambda = Node::toLmDType(_max_lambda);

      // First, solve it as a reduction at this lambda
      initSetsAtMaxLambda();
      

      // Now, for each of these sets, run a 
      

    }

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
