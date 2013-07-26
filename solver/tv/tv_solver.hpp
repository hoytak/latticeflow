#ifndef _TV_SOLVERS_H_
#define _TV_SOLVERS_H_

#define _DEBUG_IN_TV_SOLVER_INCLUDES_

#include "../common.hpp"
#include "tv_push_relabel.hpp"
#include "tv_regpath_node.hpp"
#include "../lattices/kernellattice.hpp"
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
  
  // We're simply using one version of the Parametric flow node as our
  // one here 
  template <class Kernel, typename dtype> using TVFlowNode = 
    PFFlowNode<Kernel, dtype, PFScaledUnweightedNodePolicy>;

  template <typename _Kernel, typename dtype = long> class TVSolver 
  {
  public:

    typedef _Kernel Kernel;
    typedef TVFlowNode<Kernel, dtype> Node;
    typedef KernelLattice<Node, Kernel::n_dimensions, Kernel> Lattice;
    typedef TV_PRFlow<dtype, Lattice> PRSolver;

    static constexpr int n_dimensions = Lattice::n_dimensions;

    typedef typename Lattice::index_vect index_vect;
    typedef typename PRSolver::node_ptr node_ptr;    
    typedef TVRegPathSegment<dtype, PRSolver> _TVRegPathSegment;

    typedef typename _TVRegPathSegment::Mode RPSMode;

    typedef LatticeArray<double, n_dimensions + 1> FuncPathArray;
    typedef LatticeArray<double, n_dimensions> FuncArray;

    static_assert(Kernel::is_geocut_applicable,
                  "Kernel is not valid for GeoCuts or TV Minimization.");

    static_assert(Kernel::n_dimensions == 2, 
                  "Currently only dealing with 2d stuff.");

    static_assert(8*sizeof(dtype) >= 32, 
                  "This class does not work with lower precision data type; use dtype of 32 bits or more.");

    TVSolver(index_vect dimensions)
      : lattice(dimensions)
      , function_set(false)
      , v_min(0)
      , v_max(0)
      , solver(lattice)
      , calculated_lambda(-1)
      , max_lambda(0)
      , node_map_at_lambda_max(lattice.shape())
    {
    }

    void setup(double *function, double _max_lambda) {
      // reset();

      v_min = *min_element(function, function + lattice.sizeWithinBounds()) - 1e-32;
      v_max = *max_element(function, function + lattice.sizeWithinBounds()) + 1e-32;

      max_lambda = max(_max_lambda, 0.1);
      
      typename Node::template NodeFiller<Lattice> filler(lattice);

      for(auto node_it = lattice.vertexIterator(); !node_it.done(); ++node_it) {
        dtype uf = toFVDType(function[node_it.nodeIndex()]);
        filler.addE1(node_it.node(), 0, uf);
        node_it.node()->pullBaseFunctionValueFromLattice();
      }

      // Since we are scaling the function value, scale the regularizer as well
      double global_scale_value = 0.5 * (2.0 / (v_max - v_min)) / max_lambda;

      for(auto edge_it = lattice.edgeIterator(); !edge_it.done(); ++edge_it) {

        dtype pwf = toUnshiftedFVDType(global_scale_value * edge_it.geocutEdgeWeight());
        
        filler.addE2(edge_it.node1(), edge_it.node2(), edge_it.edgeIndex(), 0, pwf, pwf, 0);
      }
    }

    void run() {

      dtype initial_lambda = Node::toScaleDType(1);

      ParametricFlowSolver<dtype, Lattice> pfs(lattice);
      
      // for(auto node_it = lattice.vertexIterator(); !node_it.done(); ++node_it) {
      //   cout << "f(" << node_it.coords() 
      //        << "; fv = " << node_it.node()->fv()
      //        << "; cfv = " << node_it.node()->cfv() 
      //        << "; cfv_predict = " << node_it.node()->cfv_predict()
      //        << endl;
      // }

      vector<vector<node_ptr> > levelset_maps = pfs.run();

      // cout << "After running!!!! " << endl;

      // for(auto node_it = lattice.vertexIterator(); !node_it.done(); ++node_it) {
      //   cout << "f(" << node_it.coords() 
      //        << "; fv = " << node_it.node()->fv()
      //        << "; cfv = " << node_it.node()->cfv() 
      //        << "; cfv_predict = " << node_it.node()->cfv_predict()
      //        << endl;
      // }

      // Now build the whole regularization path
      _constructInitialRegPathsFromSolvedLattice(levelset_maps, initial_lambda);
    }

    void run(double *function, double max_lambda) {
      setup(function, max_lambda);
      run();
    }

    FuncArray runSingleLambda(double *function, double lambda) {

      setup(function, lambda);
      
      ParametricFlowSolver<dtype, Lattice> pfs(lattice);

      pfs.run();
      
      FuncArray R(lattice.shape());
      
      for(auto it = lattice.vertexIterator(); !it.done(); ++it) {
        R[it.coords()] = toFValue(Node::translateFromRaw(it->r()));
      }

      return R;
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Retrieving the path
    
    FuncPathArray getRegularizationPath(const vector<double>& _lambda_values) {

      size_t n_lambda = _lambda_values.size();
      
      struct LambdaIdx {
        dtype value; 
        size_t index;
        bool operator<(const LambdaIdx& lm) const { return value > lm.value; }
      };

      vector<LambdaIdx> lambda_idx_values(n_lambda);
      vector<dtype> lambda_calc_values(n_lambda);

      for(size_t i = 0; i < n_lambda; ++i) 
        lambda_idx_values[i] = {Node::toScaleDType(_lambda_values[i] / max_lambda), i};

      sort(lambda_idx_values.begin(), lambda_idx_values.end());

      for(size_t i = 0; i < n_lambda; ++i) 
        lambda_calc_values[i] = lambda_idx_values[i].value;
      
      if(calculated_lambda == -1
         || lambda_calc_values.front() > calculated_lambda 
         || lambda_calc_values.back() < 0) {
        throw out_of_range("Given lambda values not in calculated range.");
      }
      
      FuncPathArray values(concat(n_lambda, lattice.shape()));

      for(auto it = lattice.vertexIterator(); !it.done(); ++it) {

        vector<dtype> path_values = traceRegularizationPath(it.coords(), lambda_calc_values);
        
        assert_equal(path_values.size(), n_lambda);

        for(size_t i = 0; i < path_values.size(); ++i)
          values[concat(lambda_idx_values[i].index, it.coords())] 
            = toFValue(path_values[i]);
      }

      return values;
    }

  private:
    ////////////////////////////////////////////////////////////////////////////////
    // Common variables

    Lattice lattice;
    bool function_set;
    double v_min, v_max;

    inline dtype toFVDType(double x) const {
      assert_leq(x / 2, v_max);
      assert_leq(v_min, x / 2);

      double centered = (x - v_min) / (v_max - v_min);
      double v = 2*centered - 1;

      assert_leq(abs(v), 2);

      return Node::toFVDType(v); // (x - (v_min + 0.5* (v_max - v_min) ) ) * (2.0 / (v_max - v_min)) );
    }

    inline dtype toUnshiftedFVDType(double x) const {
      return Node::toFVDType(x * (2.0 / (v_max - v_min)) );
    }

    inline double toFValue(dtype x) const {
      double centered = (Node::toFValue(x) + 1) / 2.0;
      double v = v_min + centered*(v_max - v_min);

      assert_equal(x, toFVDType(v));
      return v;
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Stuff dealing with the solver

    PRSolver solver;
    dtype calculated_lambda;
    double max_lambda;
    LatticeArray<_TVRegPathSegment*, n_dimensions> node_map_at_lambda_max;

    ////////////////////////////////////////////////////////////////////////////////
    // For managing new regpath instances 

    //. Use deque, as addresses are guaranteed not to change; still
    //have O(1) lookup by index (key)
    deque<_TVRegPathSegment> _regpathsegment_hold;

    _TVRegPathSegment* getNew_TVRegPathSegment() {
      _regpathsegment_hold.emplace_back
        (uint(_regpathsegment_hold.size()), lattice, solver);

      return &_regpathsegment_hold.back();
    }

    _TVRegPathSegment* lookupRPSFromKey(uint key) {
      return &(_regpathsegment_hold[key]);
    }

    ////////////////////////////////////////////////////////////////////////////////
    // INITIAL: methods for constructing the initial reg paths

    void _constructInitialRegPathsFromSolvedLattice(
        const vector<vector<node_ptr> >& level_sets, dtype solved_lamba) {

#ifndef NDEBUG
        for(auto& n : lattice)
          assert(n.keyIsClear());
#endif

      ////////////////////////////////////////////////////////////////////////////////
      // Build an initial set of regions at the primary lambda value

      vector<_TVRegPathSegment*> initial_rps_nodes;

      for(auto& v : node_map_at_lambda_max)
        v = nullptr;

      for(const auto& nv : level_sets) {
        _TVRegPathSegment *rps = getNew_TVRegPathSegment();        

        rps->setupAsInitial(nv.begin(), nv.end(), solved_lamba);
        
        initial_rps_nodes.push_back(rps);

        for(node_ptr n : nv) {
          n->template setKey<0>(rps->ci().key);
          assert(rps == lookupRPSFromKey(n->key()));

          node_map_at_lambda_max[lattice.getCoords(n)] = rps;
        }
      }

#ifndef NDEBUG
        for(auto rps : node_map_at_lambda_max) 
          assert(rps != nullptr);
#endif

      ////////////////////////////////////////////////////////////////////////////////
      // Now go through and construct all the neighborhood maps.  
      auto register_pair = [&](node_ptr n1, node_ptr n2) {
        _TVRegPathSegment *rps1 = node_map_at_lambda_max[lattice.getCoords(n1)];
        _TVRegPathSegment *rps2 = node_map_at_lambda_max[lattice.getCoords(n2)];
        
        rps1->neighbors().insert(rps2);
        rps2->neighbors().insert(rps1);
      };
      
      solver.constructNeighborhoodKeyPairs(register_pair);
      

      ////////////////////////////////////////////////////////////////////////////////
      // Now build the entire path

      struct FunPoint {
        dtype lambda;
        enum {Join, Split, SplitUB} mode;
        _TVRegPathSegment *rps1, *rps2;
        
        bool operator<(const FunPoint& jp) const {return lambda < jp.lambda;}
      }; 

      priority_queue<FunPoint> run_heap;

      ////////////////////////////////////////////////////////////////////////////////
      // Convenience functions to enable quick registration of segments

      auto registerPossibleJoin = 
        [&run_heap](_TVRegPathSegment *rps1, _TVRegPathSegment *rps2, dtype lambda_start) {

        cout << "REGISTERING JOIN" << endl;

        dtype join_lambda = _TVRegPathSegment::calculateJoins(rps1, rps2, lambda_start);
        if(join_lambda > 0)
          run_heap.push(FunPoint({join_lambda, FunPoint::Join, rps1, rps2}));
      };

      auto registerPossibleSplit = 
        [&run_heap](_TVRegPathSegment *rps, dtype lambda_start) {

        auto sp_info = rps->calculateSplit(lambda_start);
        
        if(sp_info.split_occurs) {
          run_heap.push(FunPoint({sp_info.split_lambda, FunPoint::Split, rps, nullptr}));
        } else if(sp_info.split_ub != 0) {
          run_heap.push(FunPoint({sp_info.split_ub, FunPoint::SplitUB, rps, nullptr}));
        }
      };

      ////////////////////////////////////////////////////////////////////////////////
      // INIT the paths for lookup

      // Init the priority queue for the joins.
      for(_TVRegPathSegment* rps1 : initial_rps_nodes) {
        for(_TVRegPathSegment* rps2 : rps1->neighbors()) {
          if(rps1 < rps2) 
            registerPossibleJoin(rps1, rps2, solved_lamba);
        }
        
        // Add them to the active lines
      }

      // Now go through and calculate all the splits, putting the
      // split upper bounds and the splits on the heap
      for(_TVRegPathSegment* rps : initial_rps_nodes) 
        registerPossibleSplit(rps, solved_lamba);

      // They are all in the priority heap at this point
      initial_rps_nodes.clear();

      ////////////////////////////////////////////////////////////////////////////////
      // Now build the reg path
                                                        
      auto activateNewLine = 
        [&, run_heap](_TVRegPathSegment* rps, dtype lambda_start) {
        
        assert_equal(lambda_start, rps->rhs_lambda);

        // Check out the joins
        for(_TVRegPathSegment* rps2 : rps->neighbors()) 
          registerPossibleJoin(rps, rps2, lambda_start);

        // Calculate potential splits
        registerPossibleSplit(rps, lambda_start);
      };

      auto isStillValid = [&](_TVRegPathSegment *rps, dtype lambda) {
        assert(rps != nullptr);
        if(rps->lhs_lambda != -1) {
          assert_leq(lambda, rps->lhs_lambda);
          assert(rps->lhs_mode != _TVRegPathSegment::Unset);
          return false;
        } else {
          return true;
        }
      };

      // Now we have something for it
      while(!run_heap.empty()) {
        FunPoint fp = run_heap.top();
        run_heap.pop();

        dtype current_lambda = fp.lambda;

        if(!isStillValid(fp.rps1, current_lambda) 
           || (fp.rps2 != nullptr && !isStillValid(fp.rps2, current_lambda)))
          continue;

        switch(fp.mode) {
        case FunPoint::Join: 
          {
            _TVRegPathSegment* rps1 = fp.rps1;
            _TVRegPathSegment* rps2 = fp.rps2;

            _TVRegPathSegment* new_rps = getNew_TVRegPathSegment();

            new_rps->setupFromJoin(rps1, rps2, current_lambda);

            rps1->deactivate();
            rps2->deactivate();

            activateNewLine(new_rps, current_lambda);

            break;
          }
        case FunPoint::Split:
          {
            _TVRegPathSegment* rps = fp.rps1;
            assert(fp.rps2 == nullptr);

            assert_equal(current_lambda, rps->ci().lambda_of_split);

            _TVRegPathSegment* new_rps1 = getNew_TVRegPathSegment();
            _TVRegPathSegment* new_rps2 = getNew_TVRegPathSegment();

            rps->applySplit(new_rps1, new_rps2, current_lambda, 
                            [&, this](uint k) {return this->lookupRPSFromKey(k); });

            rps->deactivate();

            activateNewLine(new_rps1, current_lambda);
            activateNewLine(new_rps2, current_lambda);

            break;
          }
        case FunPoint::SplitUB:
          {
            _TVRegPathSegment* rps = fp.rps1;          
            assert(fp.rps2 == nullptr);

            assert_equal(current_lambda, rps->ci().split_calculation_done_to_lambda);

            registerPossibleSplit(rps, current_lambda);

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
        assert_equal(_regpathsegment_hold.back().lhs_mode, _TVRegPathSegment::Unset);

        for(const auto& rps : _regpathsegment_hold) {
          if(&rps != &(_regpathsegment_hold.back())) {
            assert_gt(rps.lhs_lambda, 0);
            assert(rps.lhs_mode != _TVRegPathSegment::Unset);
          }
        }
      }

      auto& last_rps = _regpathsegment_hold.back();
      last_rps.lhs_lambda = 0;
      last_rps.lhs_mode = _TVRegPathSegment::Initial;
      last_rps.deactivate();

      // And we are done.
      calculated_lambda = solved_lamba;
    }

    vector<dtype> traceRegularizationPath(const index_vect& pos, const vector<dtype>& lambdas) const {
      node_ptr n = lattice(pos);

      assert(lattice.withinBounds(pos));
      assert(node_map_at_lambda_max.withinBounds(pos));

      _TVRegPathSegment* cur_rps = node_map_at_lambda_max[pos];

      if(cur_rps == nullptr) {
        assert(!lattice.withinBounds(pos));
        return vector<dtype>(lambdas.size(), 0);
      }

      size_t idx = 0;
      vector<dtype> r_values(lambdas.size());

      if(unlikely(lambdas.size() == 0))
        return r_values;

      assert_geq(cur_rps->rhs_lambda, lambdas[0]);

      // bool print = (pos == index_vect{0,0});

      // if(print)
      //   cout << "Tracing position " << pos << ":";

      while(true) {
        // cout << "P0: cur_rps->lhs_lambda = " << cur_rps->lhs_lambda
        //      << "; lambdas[idx] = " << lambdas[idx] << endl;

        while(cur_rps->lhs_lambda > lambdas[idx]) {
          // if(print)
          //   cout << "P1: cur_rps->lhs_lambda = " << cur_rps->lhs_lambda
          //        << "; lambdas[idx] = " << lambdas[idx] << endl;

          if(cur_rps->lhs_mode == _TVRegPathSegment::Split) {
            _TVRegPathSegment *rps0 = cur_rps->lhs_nodes[0];
            _TVRegPathSegment *rps1 = cur_rps->lhs_nodes[1];

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

            cur_rps = cur_rps->lhs_nodes[0];
          }
        }

        while(cur_rps->lhs_lambda <= lambdas[idx]) {
          assert_leq(cur_rps->lhs_lambda, lambdas[idx]);
          assert_leq(lambdas[idx], cur_rps->rhs_lambda);

          r_values[idx] = cur_rps->get_r_AtLambda(lambdas[idx]);

          // if(print)
          //   cout << "P2: cur_rps->lhs_lambda = " << cur_rps->lhs_lambda
          //        << "; lambdas[idx] = " << lambdas[idx]
          //        << "; r = " << r_values[idx] << endl; 
          
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

  typedef shared_ptr<LatticeArray<double, 2> > FuncMapPtr;

  template <typename Kernel, typename dtype = long>
  FuncMapPtr calculate2dTV(size_t nx, size_t ny, 
                           double *function, double lambda) {
    TVSolver<Kernel, dtype> solver({nx, ny});
    return FuncMapPtr(new LatticeArray<double, 2>(solver.runSingleLambda(function, lambda)));
  }

  typedef shared_ptr<LatticeArray<double, 3> > RegPathPtr;

  template <typename Kernel, typename dtype = long>
  RegPathPtr calculate2dTV(size_t nx, size_t ny, 
                           double *function, 
                           vector<double> lambda) {

    if(lambda.empty()) {
      cerr << "Must supply list of regularization values at which to calculate the path." << endl;
      return RegPathPtr(new LatticeArray<double, 3>({0,nx,ny}));
    }

    for(double& lm : lambda) {
      if(lm < 0) {
        cerr << "Ignoring reg value " << lm << " that is less than 0." << endl;
        lm = 0;
      }
    }

    TVSolver<Kernel, dtype> solver({nx, ny});

    solver.run(function, *max_element(lambda.begin(), lambda.end()));

    return RegPathPtr(new LatticeArray<double, 3>(solver.getRegularizationPath(lambda)));
  }

}; 

#include "../common/debug_flymake_test.hpp"

#endif /* _TV_REGPATH_H_ */






