#ifndef _TV_REGPATH_H_
#define _TV_REGPATH_H_

#include "../common.hpp"
#include "../lattices/kernellattice.hpp"
#include "../parametricflow/mn_netflow_solver.hpp"
#include "../energy_minimization/energy.hpp"
#include "tv_flow_node.hpp"
#include "tv_push_relabel.hpp"
#include "tv_regpath_node.hpp"

#include <set>
#include <list>
#include <vector>
#include <deque>
#include <queue>
#include <unordered_map>
#include <stdexcept>

namespace latticeQBP {
  
  template <typename _Kernel, typename dtype = long> class TVSolver 
  {
  public:

    typedef _Kernel Kernel;
    typedef TVFlowNode<Kernel, dtype> Node;
    typedef KernelLattice<Node, Kernel::n_dimensions, Kernel> Lattice;
    typedef TV_PRFlow<dtype, Lattice> PRSolver;

    typedef typename Lattice::index_vect index_vect;
    typedef typename PRSolver::node_ptr node_ptr;    
    typedef TVRegPathSegment<dtype, PRSolver> TVRegPathSegment;

    typedef typename TVRegPathSegment::Mode RPSMode;

    static constexpr int n_dimensions = Lattice::n_dimensions;

    static_assert(Kernel::is_geocut_applicable,
                  "Kernel is not valid for GeoCuts or TV Minimization.");

    static_assert(Kernel::n_dimensions == 2, 
                  "Currently only dealing with 2d stuff.");

    static_assert(8*sizeof(dtype) >= 32, 
                  "This class does not work with lower precision data type; use dtype of 32 bits or more.");

    TVSolver(index_vect dimensions)
      : lattice(dimensions)
      , solver(lattice)
      , calculated_lambda(-1)
      , node_map_at_lambda_max(lattice.shape())
      , function_set(false)
      , v_min(0)
      , v_max(0)
    {
    }

    void set(double *function) {
      v_min = *min_element(function, function + lattice.sizeWithinBounds());
      v_max = *max_element(function, function + lattice.sizeWithinBounds());

      for(auto it = lattice.indexIterator(); !it.done(); ++it) {
        const index_vect& idx = it.coords();
        lattice(idx)->setBaseFunctionValue(function[it.boundedIndex()]);
      }

      typename Node::template NodeFiller<Lattice> filler(lattice);

      for(auto edge_it = lattice.edgeIterator(); !edge_it.done(); ++edge_it) {
        
        index_vect c1 = edge_it.latticeCoordOf1();
        size_t idx_src = edge_it.nodeIndexOf1();
        size_t idx_dest = edge_it.nodeIndexOf2();

        dtype pwf = Node::toFVDType(0.5*edge_it.geocutEdgeWeight() 
                                    * abs(function[idx_src] - function[idx_dest]));
        
        filler.addE2(edge_it.node1(), edge_it.node2(), edge_it.edgeIndex(), 0, pwf, pwf, 0);
      }
    }

    void buildFullPathFromSolvedLambda(double initial_lambda) {

      // We assume that the lattice has been solved with each level
      // being a different partition.
      if(initial_lambda < calculated_lambda)
        return;

      // Solve the initial path; this can be swapped out later with a more efficient routine.
      for(auto& n : lattice)
        n.setFunctionValue(lattice, 0, initial_lambda);
      
      NetflowReductionSolver<dtype, Lattice> nrs(lattice);
      nrs.run();

      // Now build the whole regularization path
      _constructInitialRegPathsFromSolvedLattice(initial_lambda);
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Retrieving the path
    typedef LatticeArray<double, n_dimensions + 1> FuncPathArray;
    
    template <typename FunctionConversionFunction>
    FuncPathArray getRegularizationPath(const vector<double>& _lambda_values,
                                        const FunctionConversionFunction& toFunc) const {

      size_t n_lambda = _lambda_values.size();

      struct LambdaIdx {
        dtype value; 
        size_t index;
        bool operator<(const LambdaIdx& lm) const { return value > lm.value; }
      };

      vector<LambdaIdx> lambda_idx_values(n_lambda);
      vector<dtype> lambda_calc_values(n_lambda);

      for(size_t i = 0; i < n_lambda; ++i) 
        lambda_idx_values[i] = {Node::toLmDType(_lambda_values[i]), i};

      sort(lambda_idx_values.begin(), lambda_idx_values.end());

      for(size_t i = 0; i < n_lambda; ++i) 
        lambda_calc_values[i] = lambda_idx_values[i].value;
      
      if(calculated_lambda == -1
         || lambda_calc_values.front() > calculated_lambda 
         || lambda_calc_values.back() < 0) {
        throw out_of_range("Given lambda values not in calculated range.");
      }

      FuncPathArray values(concat(n_lambda, lattice.shape()));

      for(IndexIterator<n_dimensions> it(lattice.shape()); !it.done(); ++it) {

        vector<dtype> path_values = traceRegularizationPath(it.coords(), lambda_calc_values);
        
        assert_equal(path_values.size(), n_lambda);

        for(size_t i = 0; i < path_values.size(); ++i)
          values[concat(lambda_idx_values[i].index, it.coords())] = toFunc(path_values[i]);
      }
    }

  private:
    ////////////////////////////////////////////////////////////////////////////////
    // Common variables

    Lattice lattice;
    bool function_set;
    double v_min, v_max;

    inline dtype toFVDtype(double x) const {
      return Node::toFVDType( (x - (v_min + 0.5* (v_max - v_min) ) ) * (2.0 / (v_max - v_min)) );
    }

    inline double toFValue(dtype x) const {
      return (0.5*(Node::toFValue(x) + 1))*(v_max - v_min) + v_min;
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Stuff dealing with the solver

    PRSolver solver;
    dtype calculated_lambda;
    LatticeArray<TVRegPathSegment*, n_dimensions> node_map_at_lambda_max;

    ////////////////////////////////////////////////////////////////////////////////
    // For managing new regpath instances 

    //. Use deque, as addresses are guaranteed not to change; still
    //have O(1) lookup by index (key)
    deque<TVRegPathSegment> _regpathsegment_hold;

    TVRegPathSegment* getNewTVRegPathSegment() {
      _regpathsegment_hold.emplace_back
        (uint(_regpathsegment_hold.size()), lattice, solver);

      return &_regpathsegment_hold.back();
    }

    TVRegPathSegment* lookupRPSFromKey(uint key) {
      return &(_regpathsegment_hold[key]);
    }

    ////////////////////////////////////////////////////////////////////////////////
    // INITIAL: methods for constructing the initial reg paths

    void _constructInitialRegPathsFromSolvedLattice(dtype solved_lamba) {

      if(DEBUG_MODE) {
        for(auto& n : lattice) 
          assert(n.keyIsClear());
      }

      ////////////////////////////////////////////////////////////////////////////////
      // Build an initial set of regions at the primary lambda value

      vector<TVRegPathSegment*> initial_rps_nodes;

      // The null one is for the dead nodes on the perimeters 
      TVRegPathSegment *null_rps = getNewTVRegPathSegment();

      for(auto it = node_map_at_lambda_max.indexIterator(); !it.done(); ++it) {

        node_ptr n = lattice(it);
        TVRegPathSegment *rps;

        if(unlikely(n->keyIsClear())) {

          if(solver.nodeIsOrphan(n)) {
            rps = null_rps;
          } else {

            dtype lvl = n->level();

            vector<node_ptr> region = solver.walkConnectedRegion
              (n, [lvl](node_ptr nn) { return abs(lvl - nn->level()) <= 1; } );

            rps = getNewTVRegPathSegment();
        
            rps->setupAsInitial(region.begin(), region.end(), solved_lamba);

            assert(rps == lookupRPSFromKey(n->key()));

            initial_rps_nodes.push_back(rps);
          }

        } else {
          rps = lookupRPSFromKey(n->key());
          
          // Create the new region; setting all the sections to the
          // correct key. 
        }

        node_map_at_lambda_max[it] = rps;
      }

      ////////////////////////////////////////////////////////////////////////////////
      // Now go through and construct all the neighborhood maps.  
      auto register_pair = [&](node_ptr n1, node_ptr n2) {
        TVRegPathSegment *rps1 = node_map_at_lambda_max[n1 - lattice.begin()];
        TVRegPathSegment *rps2 = node_map_at_lambda_max[n2 - lattice.begin()];
        
        rps1->neighbors().insert(rps2);
        rps2->neighbors().insert(rps1);
      };
      
      solver.constructNeighborhoodKeyPairs(register_pair);
      

      ////////////////////////////////////////////////////////////////////////////////
      // Now build the entire path

      struct FunPoint {
        dtype lambda;
        enum {Join, Split, SplitUB} mode;
        TVRegPathSegment *rps1, *rps2;
        
        bool operator<(const FunPoint& jp) const {return lambda < jp.lambda;}
      }; 

      priority_queue<FunPoint> run_heap;

      dtype current_lambda = solved_lamba;

      ////////////////////////////////////////////////////////////////////////////////
      // Convenience functions to enable quick registration of segments

      auto registerPossibleJoin = 
        [&, current_lambda](TVRegPathSegment *rps1, TVRegPathSegment *rps2) {

        dtype join_lambda = TVRegPathSegment::calculateJoins(rps1, rps2, current_lambda);
        if(join_lambda > 0)
          run_heap.push(FunPoint({join_lambda, FunPoint::Join, rps1, rps2}));
      };

      auto registerPossibleSplit = 
        [&, current_lambda](TVRegPathSegment *rps) {

        auto sp_info = rps->calculateSplit(current_lambda);
        
        if(sp_info.split_occurs) {
          run_heap.push(FunPoint({sp_info.split_lambda, FunPoint::Split, rps, nullptr}));
        } else if(sp_info.split_ub != 0) {
          run_heap.push(FunPoint({sp_info.split_ub, FunPoint::SplitUB, rps, nullptr}));
        }
      };

      ////////////////////////////////////////////////////////////////////////////////
      // INIT the paths for lookup

      // Init the priority queue for the joins.
      for(TVRegPathSegment* rps1 : initial_rps_nodes) {
        for(TVRegPathSegment* rps2 : rps1->neighbors()) {
          if(rps1 < rps2) 
            registerPossibleJoin(rps1, rps2);
        }
        
        // Add them to the active lines
      }

      // Now go through and calculate all the splits, putting the
      // split upper bounds and the splits on the heap
      for(TVRegPathSegment* rps : initial_rps_nodes) 
        registerPossibleSplit(rps);

      // They are all in the priority heap at this point
      initial_rps_nodes.clear();

      ////////////////////////////////////////////////////////////////////////////////
      // Now build the reg path
                                                        
      auto activateNewLine = 
        [&, run_heap, current_lambda](TVRegPathSegment* rps) {
        
        assert_equal(current_lambda, rps->rhs_lambda);

        // Check out the joins
        for(TVRegPathSegment* rps2 : rps->neighbors()) 
          registerPossibleJoin(rps, rps2);

        // Calculate potential splits
        registerPossibleSplit(rps);
      };

      auto isStillValid = [&, current_lambda](TVRegPathSegment *rps) {
        assert(rps != nullptr);
        if(rps->lhs_lambda != -1) {
          assert_geq(current_lambda, rps->lhs_lambda);
          assert(rps->lhs_mode != TVRegPathSegment::Unset);
          return false;
        } else {
          return true;
        }
      };

      // Now we have something for it
      while(!run_heap.empty()) {
        FunPoint fp = run_heap.top();
        run_heap.pop();

        assert_leq(fp.lambda, current_lambda);
        current_lambda = fp.lambda;

        if(!isStillValid(fp.rps1) || (fp.rps2 != nullptr && !isStillValid(fp.rps2)))
          continue;

        switch(fp.mode) {
        case FunPoint::Join: 
          {
            TVRegPathSegment* rps1 = fp.rps1;
            TVRegPathSegment* rps2 = fp.rps2;

            TVRegPathSegment* new_rps = getNewTVRegPathSegment();

            new_rps->setupFromJoin(current_lambda, rps1, rps2);

            rps1->deactivate();
            rps2->deactivate();

            activateNewLine(new_rps);

            break;
          }
        case FunPoint::Split:
          {
            TVRegPathSegment* rps = fp.rps1;
            assert(fp.rps2 == nullptr);

            assert_equal(current_lambda, rps->constructionInfo()->lambda_of_split);

            TVRegPathSegment* new_rps1 = getNewTVRegPathSegment();
            TVRegPathSegment* new_rps2 = getNewTVRegPathSegment();

            rps->applySplit(new_rps1, new_rps2, current_lambda, 
                            [&, this](uint k) {return lookupRPSFromKey(k); });

            rps->deactivate();

            activateNewLine(new_rps1);
            activateNewLine(new_rps2);

            break;
          }
        case FunPoint::SplitUB:
          {
            TVRegPathSegment* rps = fp.rps1;          
            assert(fp.rps2 == nullptr);

            assert_equal(current_lambda, 
                         rps->constructionInfo()->split_calculation_done_to_lambda);

            registerPossibleSplit(rps);

            break;
          }
        }
      }

      // At this point, everything should be solved.  If we are
      // correct, there should have been a join BEFORE lambda = 0, as
      // the lambda term is rounded up in the calculations. Thus we
      // should have exactly that: the last RegPath added has no
      // lhs_lambda value, its rhs_lambda value is > 0, and all the
      // other RegPath entries have lhs_lambda > 0. 

      if(DEBUG_MODE) {
        assert_equal(_regpathsegment_hold.back().lhs_lambda, -1);
        assert_equal(_regpathsegment_hold.back().lhs_mode, TVRegPathSegment::Unset);

        for(const auto& rps : _regpathsegment_hold) {
          if(&rps != &(_regpathsegment_hold.back())) {
            assert_gt(rps.lhs_lambda, 0);
            assert(rps.lhs_mode != TVRegPathSegment::Unset);
          }
        }
      }

      auto& last_rps = _regpathsegment_hold.back();
      last_rps.lhs_lambda = 0;
      last_rps.lhs_mode = TVRegPathSegment::Initial;
      last_rps.deactivate();

      // And we are done.
      calculated_lambda = solved_lamba;
    }

    vector<dtype> traceRegularizationPath(const index_vect& pos, const vector<dtype>& lambdas) const {
      node_ptr n = lattice(pos);

      TVRegPathSegment* cur_rps = node_map_at_lambda_max[pos];

      size_t idx = 0;
      vector<dtype> r_values(lambdas.size());

      if(unlikely(lambdas.size() == 0))
        return r_values;

      assert_geq(cur_rps->rhs_lambda, lambdas[0]);

      while(true) {
        while(cur_rps->lhs_lambda > lambdas[idx]) {
          if(cur_rps->lhs_mode == TVRegPathSegment::Split) {
            TVRegPathSegment *rps0 = cur_rps->lhs_nodes[0];
            TVRegPathSegment *rps1 = cur_rps->lhs_nodes[1];

            assert(rps0 != nullptr);
            assert(rps1 != nullptr);

            if(rps0->containsNode(n)) {
              cur_rps = rps0;
            } else {
              assert(rps1->containsNode(n));
              cur_rps = rps1;
            }

          } else {
            assert(cur_rps->lhs_nodes[0] != nullptr);
            assert(cur_rps->lhs_nodes[1] == nullptr);

            cur_rps->lhs_nodes[0];
          }
        }

        while(cur_rps->lhs_lambda <= lambdas[idx]) {
          r_values[idx] = cur_rps->get_r_AtLambda(lambdas[idx]);
          ++idx;
          
          if(idx == lambdas.size())
            goto r_values_finished;
        }
      }
    r_values_finished:;
      return r_values;
    }
  };

  ////////////////////////////////////////////////////////////////////////////////
  // A one-shot function to reference against

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

#ifdef EMACS_FLYMAKE

#include "../kernels/kernels.hpp"

namespace latticeQBP {
  template class TVSolver<Full2d_4, long>;
};
#endif


#endif /* _TV_REGPATH_H_ */
