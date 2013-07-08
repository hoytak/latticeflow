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
    typedef typename CompType<dtype>::Type comp_type;

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

    typedef PRFlow<dtype, Lattice, 0> PRSolver;

    PRSolver pr_solver;

    ////////////////////////////////////////////////////////////////////////////////
    // Initializing all the sets 

    struct NodePivot;
    typedef shared_ptr<set<node_ptr> > nodeset_ptr;

    struct RegPathSegment {
      typedef enum  {Terminal, Join, Split} Mode ;

      Mode rhs_mode, lhs_mode; 

      dtype rhs_lambda, lhs_lambda;
      dtype rhs_r;
      dtype r_at_zero;

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
    
    dtype pathIntersection(RegPathSegment* rhs1, RegPathSegment* rhs2) const {
      // This function returns the lambda value at which the two paths
      // join, or zero if they never cross.
      return 0;
    }

    ////////////////////////////////////////////////////////////
    // Adding nodes in there 

    typedef typename RegPathSegment::Mode RPSMode;

    template <typename ForwardIterator> 
    void addPivot(ForwardIterator it_start, ForwardIterator it_end, dtype node_rhs_lambda,
                  RPSMode rhs_mode)
    {
      
      RegPathSegment* rps = getNewRegPathSegment();

      // First, go through and check all the other paths to see if any
      // will join this one.  If so, then 





      ////////////////////////////////////////////////////////////////////////////////
      // First, see if it's possible to solve the transhipment problem
      // at the lowest possible lambda

      typedef typename PRSolver::partitioninfo_ptr partitioninfo_ptr;

      dtype min_lambda = known_join_lambda;
      vector<partitioninfo_ptr> current_partition_info_split;

      while(true) { 

        StableAverage<dtype> fv_avg;
      
        for(ForwardIterator it = it_start; it != it_end; ++it) {
          node_ptr n = (*it);
          n->setFunctionValue(lattice, 0, min_lambda); 
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
        vector<partitioninfo_ptr> piv = pr_solver.cleanupAndGetBlocks(it_start, it_end, key);      
        // See if we've found the correct lower bound on the lambda
        if(piv.size() == 1)
          break;
        
        // Okay, now recallibrate the lambda needed to correspond to the minimum breaking point. 
        current_partition_info_split = piv;

        // Need a set of equations to calculate the 

        struct PartInfo {
          partitioninfo_ptr pt;
          dtype qii_sum, gamma_sum;
        };

        // calculate the slope and intercept terms 
        size_t R_size = 0;
        for(const partitioninfo_ptr& pt : piv)
          R_size += pt->nodes.size();

        vector<PartInfo> p_info(piv.size());

        for(size_t i = 0; i < piv.size(); ++i) {
          const partitioninfo_ptr& pt = piv[i];
          PartInfo& pi = p_info[i];

          pi = {pt, 0, 0};

          dtype lm_qii_sum = 0, qii_sum = 0, total_val = 0;

          for(node_ptr n : piv->nodes) {

            // Subtracting by fv_avg is for numerical stability
            lm_qii_sum += n->cfv() - fv_avg;

            pi.qii_sum += n->fv();

            total_val += n->fv_predict() - fv_avg;
          }
          
          pi.gamma_sum = (total_val
                          - lm_qii_sum
                          - (piv->is_on ? 1 : -1)*piv->cut_value);

        }

        // Get the rest of the components to calculate the shape
        comp_type qii_total = 0;
        comp_type gamma_total = 0;

        for(const PartInfo& pi : p_info) {
          qii_total += pi.qii_sum;
          gamma_total += pi.gamma_sum;
        }

        // Now go through and see which one has the largest lambda 
        dtype max_lambda_so_far = 0;

        for(const PartInfo& pi : p_info) {
          comp_type lambda_coeff = R_size * comp_type(pi.qii_sum)   - pi.size * qii_total;
          comp_type lambda_intcp = R_size * comp_type(pi.gamma_sum) - pi.size * gamma_total;
          comp_type cut = R_size * comp_type(pi.pt->cut_value);

          dtype calc_lambda; 

          // is_on being true means that this partition was on the
          // high end, so at the lambda = 0 end, it's got too much
          // flow if this is the blocking cut section.  This means
          // that the incoming flow must decrease with increasing
          // lambda, and that the original intercept term must be
          // positive.  Thus we are looking for the point where it
          // hits the cut.


          if(  (pi.pt->is_on  && ( lambda_coeff >= 0 || cut >= lambda_intcp)  )
               || (!pi.pt->is_on && ( lambda_coeff <= 0 || cut >= -lambda_intcp) ) ) {

            // This means it is not the part that contains the
            // blocking flow. 
            continue;
          }
            
          calc_lambda = Node::getLambdaFromQuotient(abs(lambda_intcp) - cut, abs(lambda_coeff));

          assert_leq(calc_lambda, node_rhs_lambda);

          if(DEBUG_MODE && piv.size() == 2 && &pi == &(p_info[1])) {
            // These should be approximately the same 
            assert_equal(calc_lambda, max_lambda_so_far);
          }
          
          max_lambda_so_far = max(calc_lambda, max_lambda_so_far);
        }
        
        // Now, attempt to recalculate 
        min_lambda = max_lambda_so_far + 1;

      } // Go back and try the cut again
      
      // At this point, we know that the 
      

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
