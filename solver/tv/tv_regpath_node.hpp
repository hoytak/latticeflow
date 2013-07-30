#ifndef _TV_REGPATH_NODE_H_
#define _TV_REGPATH_NODE_H_

#include "../common.hpp"
#include "tv_push_relabel.hpp"

#include <set>
#include <map>
#include <vector>

#define PRINT_SPLIT_MESSAGES false
#define PRINT_PATH_MESSAGES false
#define PRINT_INTERATION_INFO_MESSAGES false
#define ENHANCE_ACCURACY false
#define PRINT_SPLIT_STAT_CHARS true

// #define PRINT_SPLIT_MESSAGES true
// #define PRINT_PATH_MESSAGES true
// #define PRINT_INTERATION_INFO_MESSAGES true

namespace latticeQBP {

// #ifdef NDEBUG
// #undef NDEBUG
// #endif

// #include "../common/debug.hpp"


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
  // Note that it is possible for lhs_lambda == rhs_lambda.  If this
  // happens, then the split / join is of more than one set at the
  // same lambda.  This is a rare occurance, so don't optimize for it.
  //
  // While things are building, the nodeset contains the current
  // nodes for both the join and split cases.  But it's then
  // cleared out when the lhs is decided. 

  
  using namespace std;

  template <typename dtype, typename TV_PR_Class> class TVRegPathSegment {
  public:
    typedef typename TV_PR_Class::Lattice Lattice;
    typedef typename TV_PR_Class::Node Node;    
    typedef typename TV_PR_Class::node_ptr node_ptr;    
    typedef typename CompType<dtype>::Type comp_type;

    typedef enum  {Unset, Initial, Join, Split} Mode;    

    TVRegPathSegment(uint _key, Lattice& _lattice, TV_PR_Class& _solver)
      : rhs_mode(Unset)
      , lhs_mode(Unset)
      , rhs_lambda(0)
      , lhs_lambda(-1)
      , adjusted_r_at_1(0)
      , adjusted_r_at_0(0)
      , _construction_info(new ConstructionInfo(_key, _lattice, _solver))
      , n_nodes(0)
      , rhs_nodes(nullptr)
      , lhs_nodes(nullptr)
    {}

    template <typename ForwardIterator>
    void setupAsInitial(const ForwardIterator& start, 
                        const ForwardIterator& end,
                        dtype lambda) {
      assert(nodeset.empty());

      assert(ci().nodeset.empty());
      ci().nodeset = vector<node_ptr>(start, end);
      assert(!ci().nodeset.empty());

      if(ENABLE_EXPENSIVE_CHECKS)
        assert_equal(set<node_ptr>(start, end).size(), ci().nodeset.size());

      rhs_lambda   = lambda;
      rhs_mode     = Initial;
      rhs_nodes[0] = nullptr;
      rhs_nodes[1] = nullptr;

      syncInformation(lambda);
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Internal data 

    ////////////////////////////////////////////////////////////////////////////////
    // Stuff for tracking things 
    Mode rhs_mode, lhs_mode;
    dtype rhs_lambda, lhs_lambda;

    ////////////////////////////////////////////////////////////////////////////////
    // DIRECTIONAL INFORMATION for testing 

  private: 

    // These are shifted by lambda_precision_bits to maintain some
    // additional numerical accuracy.
    dtype adjusted_r_at_1;
    dtype adjusted_r_at_0;

    static constexpr int n_bits_adjust = Node::n_bits_function_room;

    template <typename T>
    static inline T adjust_r(T r) { 
      const int _n_bits_adjust = n_bits_adjust;
      return r * (T(1) << _n_bits_adjust); 
    }

    template <typename T>
    static inline T deadjust_r(T r) { 
      const int _n_bits_adjust = n_bits_adjust;
      
      r += (r < 0 ? -1 : 1) * (T(1) << (_n_bits_adjust - 1));

      return r / (T(1) << _n_bits_adjust);
    }

  public:
    inline dtype get_r_AtLambda(dtype lambda) const {

      return deadjust_r(adjusted_r_at_0
                        + Node::multFVScale(adjusted_r_at_1 - adjusted_r_at_0, lambda));
    }

  private:
    struct RegionInformation {
      
      // These are referenced from the base zero_reference.  
      dtype lm_qii_sum_d,  r_sum_d;
      dtype zero_reference;
      dtype partition_size;

      comp_type qii_sum;

      // These are absolute values
      comp_type gamma_sum;
    }; 

    template <typename ForwardIterator>
    inline RegionInformation getRegionInfo(const ForwardIterator& start, 
                                           const ForwardIterator& end, 
                                           dtype lambda
                                           ) const {

      assert(start != end);

      RegionInformation ri = {0,0,0,0,0,0};

      ri.zero_reference = (*start)->r();

      for(ForwardIterator it = start; it != end; ++it) {

        node_ptr n = *it;

        n->_debug_checkLevelsetMethodsNode(ci().lattice);
        assert_equal(n->currentScale(), lambda);

        assert_equal(n->cfv(), n->fv(lambda));
        
        // Subtracting by fv_avg is for numerical stability
        ri.lm_qii_sum_d += n->lm_qii() - ri.zero_reference;
        ri.r_sum_d += n->r() - ri.zero_reference;

        ri.qii_sum += n->qii();

        ++ri.partition_size;
      }

      ri.gamma_sum = (ri.r_sum_d - ri.lm_qii_sum_d);
      return ri;
    }

#ifndef NDEBUG
    void checkKeySynced() const {
      ci().solver.checkPartitionedSection(ci().nodeset.begin(), ci().nodeset.end(), ci().key);
    }
#else
    void checkKeySynced() const { }
#endif

    void syncInformation(dtype lambda) {
      
      assert(!ci().nodeset.empty());

      // Update the lattice
      n_nodes = 0;
      for(node_ptr n : ci().nodeset) {
        
        n->setOffsetAndScale(ci().lattice, 0, lambda); 
        n->template setKey<0>(ci().key);
        ++n_nodes;
      }

      auto start = ci().nodeset.begin();
      auto end = ci().nodeset.end();

      RegionInformation ri = getRegionInfo(start, end, lambda);

      // Get the information from these values to get the updated
      // directional information

      // The value of 
      adjusted_r_at_0 = toDType(adjust_r(comp_type(ri.gamma_sum)) / ri.partition_size);
      adjusted_r_at_1 = toDType(adjust_r(ri.qii_sum + ri.gamma_sum) / ri.partition_size);

      // cout << "lm_qii_sum_d = " << ri.lm_qii_sum_d << endl; 
      // cout << "r_sum_d = " << ri.r_sum_d << endl;
      // cout << "zero_reference = " << ri.zero_reference<< endl;
      // cout << "partition_size = " << ri.partition_size<< endl;
      // cout << "qii_sum = " << ri.qii_sum<< endl;
      // cout << "gamma_sum = " << ri.gamma_sum<< endl;
      
      // cout << "lambda = " << lambda << endl;
      // cout << "adjusted_r_at_0 = " << adjusted_r_at_0<< endl;
      // cout << "adjusted_r_at_1 = " << adjusted_r_at_1<< endl;

      // cout << "deadjust_r(adjusted_r_at_0) = " << deadjust_r(adjusted_r_at_0) << endl;
      // cout << "deadjust_r(adjusted_r_at_1) = " << deadjust_r(adjusted_r_at_1) << endl;
      

      // Now, make sure that it's correct...

#ifndef NDEBUG

      assert_equal(deadjust_r(adjusted_r_at_0), get_r_AtLambda(Node::toScaleDType(0)));
      assert_equal(deadjust_r(adjusted_r_at_1), get_r_AtLambda(Node::toScaleDType(1)));

      dtype r_calc = get_r_AtLambda(lambda);

      dtype rhs_r = toDType( (comp_type(ri.zero_reference)*ri.partition_size 
                              + (ri.r_sum_d 
                                 + (ri.partition_size >> 2))) 
                             / ri.partition_size);

      assert_leq(abs(r_calc - rhs_r), 1);

      checkKeySynced();
#endif
    }


    ////////////////////////////////////////////////////////////////////////////////
    // CONSTRUCTION Aids

  public:
    struct ConstructionInfo {
      ConstructionInfo(uint _key, Lattice& _lattice, TV_PR_Class& _solver) 
        : key(_key)
        , lattice(_lattice)
        , solver(_solver)
        , split_calculation_done_to_lambda(-1)
        , lambda_of_split(-1)
      {}

      const uint key;
      Lattice& lattice;
      TV_PR_Class& solver;

      set<TVRegPathSegment*> neighbors;
      multimap<dtype, TVRegPathSegment*> join_points;
      vector<node_ptr> nodeset;

      // For potential split
      dtype split_calculation_done_to_lambda;

      dtype lambda_of_split;

      typename TV_PR_Class::cutinfo_ptr split_information;
    };

  private:
    ConstructionInfo* _construction_info;
    
  public:
    ConstructionInfo& ci() const {
      assert(_construction_info != nullptr);
      return *_construction_info;
    }

    set<TVRegPathSegment*>& neighbors() const {
      assert(_construction_info != nullptr);
      return _construction_info->neighbors;
    }

    void registerJoinPoint(dtype join_lambda, TVRegPathSegment* rps) const {
      ci().join_points.insert(make_pair(join_lambda, rps));
    }

    dtype firstJoinPoint() const {
      assert(!ci().join_points.empty());
      return ci().join_points.rbegin()->first;
    }

    // Conditionally activate the nodes 
    void deactivate() {

      assert(lhs_lambda != -1);

      // Depending on the rhs mode, create a set with all the active
      // nodes in it. 7
      if(rhs_mode == Split || rhs_mode == Initial) {
        nodeset.swap(ci().nodeset);
        sort(nodeset.begin(), nodeset.end());
      }

      if(DEBUG_MODE) {
        // Run a bunch of assert

        switch(rhs_mode) {
        case Split: 
          assert(rhs_nodes[0] != nullptr);
          assert(rhs_nodes[1] == nullptr);

          assert(rhs_nodes[0]->lhs_nodes[0] == this 
                 || rhs_nodes[0]->lhs_nodes[1] == this);

          assert_equal(rhs_lambda, rhs_nodes[0]->lhs_lambda);
          assert_equal(rhs_nodes[0]->lhs_mode, Split);

          break;

        case Join:
          assert(rhs_nodes[0] != nullptr);
          assert(rhs_nodes[1] != nullptr);

          assert(rhs_nodes[0]->lhs_nodes[0] == this);
          assert(rhs_nodes[0]->lhs_nodes[1] == nullptr);
          assert(rhs_nodes[0]->lhs_lambda == rhs_lambda);
          assert_equal(rhs_nodes[0]->lhs_mode, Join);

          assert(rhs_nodes[1]->lhs_nodes[0] == this);
          assert(rhs_nodes[1]->lhs_nodes[1] == nullptr);
          assert(rhs_nodes[1]->lhs_lambda == rhs_lambda);
          assert_equal(rhs_nodes[1]->lhs_mode, Join);
          break;

        case Initial:
          assert(rhs_nodes[0] == nullptr);
          assert(rhs_nodes[1] == nullptr);
          assert(rhs_lambda >= lhs_lambda);
          break;

        case Unset:
          assert(rhs_mode != Unset);
          break;
        }

        // Ensure it is not in any neighborhood maps
#ifndef NDEBUG        
        for(TVRegPathSegment* rps : ci().neighbors) {
          assert(rps->neighbors().find(this) 
                 == rps->neighbors().end());
        }
        
      if(ENABLE_EXPENSIVE_CHECKS)
        for(node_ptr n = ci().lattice.begin(); n != ci().lattice.end(); ++n)
          assert_notequal(n->key(), ci().key);
#endif
      }

      assert(_construction_info != nullptr);

#ifndef NDEBUG
      cout << "Deactivated node " << ci().key << "." << endl;
#endif

      delete _construction_info;
      _construction_info = nullptr;
    }


    ////////////////////////////////////////////////////////////////////////////////
    // SPLITS

  private:

    struct DetailedSplitInfo{ 
      bool split_occurs;
      bool lambda_is_exactly_known;

      // This gives the exact lambda of this cut; if it is one less
      // than this, then the cut will form.  Above this, there may be
      // a new cut, but at least this one won't form.
      dtype lambda_of_split_capacity;

      typename TV_PR_Class::cutinfo_ptr cut;
    };

    DetailedSplitInfo BisectionCheck(dtype lambda_lb, dtype lambda_ub) const {

      int bisect_rate = 4;

      typename TV_PR_Class::cutinfo_ptr best_cut;

      bool cut_occured = false;

      while(lambda_lb + 1 < lambda_ub) {

        dtype try_lambda = (bisect_rate != 1
                            ? ((lambda_lb + 1) + ((lambda_ub - (lambda_lb + 1)) / (dtype(1) << bisect_rate))) 
                            : (lambda_ub + lambda_lb) / 2);
        
        // Calculate the split point...
        ci().solver.setRegionToLambda(ci().nodeset.begin(), 
                                      ci().nodeset.end(), try_lambda, true);
      
        auto cut = ci().solver.runPartitionedSection(ci().nodeset.begin(), ci().nodeset.end(), ci().key);

        if(!cut->any_cut) {
          lambda_ub = try_lambda;
        } else {
          lambda_lb = try_lambda;
          best_cut = cut;
          cut_occured = true;
          if(bisect_rate > 1)
            --bisect_rate;
        }
      }

      if(cut_occured)
        return DetailedSplitInfo{true, true, lambda_lb, best_cut};
      else
        return DetailedSplitInfo{false, false, lambda_lb, nullptr};
    }

    DetailedSplitInfo _calculateSingleSplit(const dtype lambda_lb, dtype lambda_ub) const {

      // ci().solver.enableChecks();

      // cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
      if(PRINT_SPLIT_MESSAGES)
        cout << ci().key << ": Calculating split at lambda = " << lambda_lb << endl;

      // Calculate the split point...
      ci().solver.setRegionToLambda(ci().nodeset.begin(), 
                                    ci().nodeset.end(), lambda_lb, false);
      
      // Now see if it can be solved at that lambda
      auto cut = ci().solver.runPartitionedSection(ci().nodeset.begin(), ci().nodeset.end(), ci().key);
      const auto& piv = cut->partitions;

      if(!cut->any_cut) {
        // cout << "     > No split. " << endl;
        return DetailedSplitInfo({false, false, 0, nullptr});
      } else {
        assert(!piv[0]->nodes.empty());
        assert(!piv[1]->nodes.empty());
      }

      dtype s0 = cut->partitions[0]->nodes.size();
      dtype s1 = cut->partitions[1]->nodes.size();

      if(ENHANCE_ACCURACY && max(s0, s1) >= 500*min(s0, s1)) {
        cout << 'E';
        return BisectionCheck(lambda_lb, lambda_ub);
      }

      comp_type c = cut->cut_value;

      if(unlikely(c == 0)) {
        // We're dealing with two disjoint sections, so should be split immediately.
        DetailedSplitInfo dsi;

        dsi.split_occurs = true;
        dsi.lambda_is_exactly_known = true;
        dsi.cut = cut;
        dsi.lambda_of_split_capacity = lambda_ub;

        return dsi;
      }

      Array<RegionInformation, 2> p_info = {
        getRegionInfo(piv[0]->nodes.begin(), piv[0]->nodes.end(), lambda_lb),
        getRegionInfo(piv[1]->nodes.begin(), piv[1]->nodes.end(), lambda_lb),
      };

      comp_type g0 = p_info[0].gamma_sum - c;
      comp_type g1 = p_info[1].gamma_sum + c;

      comp_type q0 = p_info[0].qii_sum;
      comp_type q1 = p_info[1].qii_sum;

      // Okay, it's a lot more accurate if we flip things around so
      // that s1 is the smaller set of nodes....  For some reason...
      
      comp_type intercept = abs((s1*g0 - s0*g1) + c * (s0 + s1));
      comp_type denom     = abs(q1*s0 - q0*s1);
        
      dtype denom_pm = (s0 + s1);

      dtype lb = Node::getScaleFromQuotient_T(intercept - max(s1,s0)*(s0 + s1), denom + denom_pm);
      dtype ub = Node::getScaleFromQuotient_T(move(intercept), denom - denom_pm);

      if(unlikely(ub < lambda_lb)) {
        if(PRINT_SPLIT_STAT_CHARS)
          cout << 'U';

        // This one seems to really mess things up when it occurs
        return BisectionCheck(lambda_lb, lambda_ub);

        // if(ENHANCE_ACCURACY)
        //   return BisectionCheck(lambda_lb, lambda_ub);
        // else
        //   return DetailedSplitInfo({true, false, lambda_lb, cut});
      }

      if(unlikely(lb > lambda_ub)) {
        if(PRINT_SPLIT_STAT_CHARS)
          cout << 'L';
        if(ENHANCE_ACCURACY)
          return BisectionCheck(lambda_lb, lambda_ub);
        else
          return DetailedSplitInfo({true, false, lambda_ub, cut});
      }

      lb = max(lb, lambda_lb);
      ub = min(ub, lambda_ub);

      comp_type gs = g0 + g1;

      auto nl0 = &(cut->partitions[0]->nodes);
      auto nl1 = &(cut->partitions[1]->nodes);

      auto value = [&,g0,g1,gs,c,s0,s1,nl0,nl1,q0,q1](dtype lm) {
        comp_type lmq0 = 0, lmq1 = 0;

        if(s0 > 32)
          lmq0 = Node::multFVScale(q0, lm);
        else
          for(node_ptr n : *nl0) 
            lmq0 += n->qii(lm);

        if(s1 > 32)
          lmq1 = Node::multFVScale(q1, lm);
        else
          for(node_ptr n : *nl1) 
            lmq1 += n->qii(lm);

        auto fl_avg = floorAverage_T(gs + lmq0 + lmq1, s0 + s1);

        return min((lmq1 + g1 - c - s1 * fl_avg), (lmq0 + g0 + c - s0 * fl_avg));
      };
                                              
      comp_type value_lb = value(lb);
      comp_type value_ub = value(ub);
                                              
      if(unlikely(value_lb >= 0)) {
        if(PRINT_SPLIT_STAT_CHARS)
          cout << 'v';
        if(ENHANCE_ACCURACY)
          return BisectionCheck(lambda_lb, lambda_ub);
        else
          return DetailedSplitInfo({true, false, lb, cut});
      }

      if(unlikely(value_ub < 0)) { 
        if(PRINT_SPLIT_STAT_CHARS)
          cout << '^';
        if(ENHANCE_ACCURACY)
          return BisectionCheck(lambda_lb, lambda_ub);
        else
          return DetailedSplitInfo({true, false, ub, cut});
      }

      while(lb + 1 < ub) {

        if(value_ub == 0)
          break;

        // cout << "Finding Value! [" << lb << ',' << ub << "] = (" << value(lb) << ',' << value(ub) << ")" << endl;
        assert_leq(value_lb, 0);
        assert_geq(value_ub, 0);
        
        dtype mid_point = (lb + ub) / 2;
        
        // double query_lm_dbl = (weight_linear * (lb + ((double(-value_lb) * (ub - lb)) / (value_ub - value_lb))) 
        //                        + (1 - weight_linear)*mid_point);

        // weight_linear *= 0.75;

        // assert_lt(query_lm_dbl, ub);
        // assert_gt(query_lm_dbl, lb); 

        dtype query_lm = mid_point; // dtype( (query_lm_dbl < mid_point) ? ceil(query_lm_dbl) : floor(query_lm_dbl));
        comp_type v_mid = value(query_lm);

        if(v_mid < 0) {
          lb = query_lm;
          value_lb = v_mid;
        } else {
          ub = query_lm;
          value_ub = v_mid;
        }
      }

      dtype pred_lambda = ub;

      // //
      // findExactPointOfCutFeasibility(max(min(calc_lambda, lambda_ub), lambda_lb), lambda_lb, cut);
      // dtype pred_lambda = min(calc_lambda, lambda_ub); 

      // if(pred_lambda < lambda_lb) {
      //   // Time to find the correct one here.  
      //   DetailedSplitInfo dsi = BisectionCheck(lambda_lb, lambda_ub);

      //   cout << "Pred too low; Correct split location = " << dsi.lambda_of_split_capacity
      //        << "; calculation = " << pred_lambda << endl; 

      //   return dsi;
      // }

      assert_geq(pred_lambda, lambda_lb);
      assert_leq(pred_lambda, lambda_ub);
        
      // cout << "DIFF: " << (pred_lambda - BisectionCheck(lambda_lb, lambda_ub)) << endl; 

      if(PRINT_SPLIT_MESSAGES)
        cout << "   " << ci().key << " Split occurs at lambda = " << pred_lambda << endl;

      return DetailedSplitInfo({true, false, pred_lambda, cut});
    }



  public:


    struct SplitInfo {
      bool split_occurs;
      dtype split_lambda;
      dtype split_ub;
    };

    SplitInfo calculateSplit(dtype current_lambda) const {
      checkKeySynced();

      if(n_nodes == 1) {
        ci().split_calculation_done_to_lambda = 0;
        ci().lambda_of_split = -1;
        
        return SplitInfo({false, -1, 0});
      }

      if(DEBUG_MODE && ci().split_calculation_done_to_lambda != -1)
        assert_leq(current_lambda, ci().split_calculation_done_to_lambda);

      // Now, go through and calculate a likely lower bound on the
      // split point.  There are situations where this ends up going
      // off, and will require recalculations, but those are few and
      // far between.
      
      dtype lambda_calc_lb = 0;

      for(auto it = ci().join_points.begin(); it != ci().join_points.end();) {
        // Remove things that are meaningless...
        dtype join_lm = it->first;
        TVRegPathSegment* rps = it->second;

        if(join_lm >= current_lambda || rps->lhs_mode != Unset) {
          auto rem_it = it;
          ++it;
          ci().join_points.erase(rem_it);
          continue;
        } else {
          if(join_lm == rps->firstJoinPoint())
            lambda_calc_lb = max(join_lm, lambda_calc_lb);

          ++it;
        }
      }
      
      DetailedSplitInfo dsi{false,false, 0}, last_dsi{false, false, 0};
      dtype lambda_calc = lambda_calc_lb;
      bool cut_exists = false;

      while(true) {
        dsi = _calculateSingleSplit(lambda_calc, current_lambda);

        if(!dsi.split_occurs)
          break;
        else 
          cut_exists = true;

        last_dsi = dsi;

        if(dsi.lambda_is_exactly_known) {
          lambda_calc = dsi.lambda_of_split_capacity;
          break;
        } else { 
          assert_gt(dsi.lambda_of_split_capacity, lambda_calc);
          assert_leq(dsi.lambda_of_split_capacity, current_lambda);

          lambda_calc = dsi.lambda_of_split_capacity + (dtype(1) << (max(0, Node::n_bits_scale_precision - 24)));

          if(lambda_calc >= current_lambda) {
            lambda_calc = current_lambda;
            break;
          }
        }
      }

      // Store that information in the computation structure
      if(cut_exists) {
        assert(last_dsi.cut != nullptr);

        assert(last_dsi.split_occurs);
        
        ci().split_information = last_dsi.cut;
        ci().lambda_of_split = last_dsi.lambda_of_split_capacity;

        // dtype true_lm = BisectionCheck(0, current_lambda).lambda_of_split_capacity;

        // cout << "Split calculated at lambda = " << ci().lambda_of_split
        //      << "; true = " << true_lm 
        //      << "; diff = " << (true_lm - ci().lambda_of_split) << endl;

        if(ENABLE_EXPENSIVE_CHECKS)
          assert_close( (1.0 - double(lambda_calc) / BisectionCheck(0, current_lambda).lambda_of_split_capacity), 0, 1e-4);

        return SplitInfo({true, lambda_calc, lambda_calc_lb});

      } else {
        ci().split_calculation_done_to_lambda = lambda_calc_lb;
        ci().lambda_of_split = -1;
        return SplitInfo({false, -1, lambda_calc_lb});
      }
    }


    template <typename RPSLookupFunction>
    void applySplit(Array<TVRegPathSegment*, 2> dest,
                    dtype lambda, 
                    const RPSLookupFunction& rpsLookup) {
      
      assert_equal(lambda, ci().lambda_of_split);
      
      // First, ensure that the cut is applied to the section so the
      // appropriate edges are saturated.
      ci().solver.applyPartioningCut(ci().split_information, ci().key);
      
      checkKeySynced();

      // Split up the nodes

      for(int i = 0; i < 2; ++i) {
        assert(dest[i]->ci().neighbors.empty());
        assert(dest[i]->ci().nodeset.empty());

        dest[i]->ci().nodeset = move(ci().split_information->partitions[i]->nodes);
        // [i]->sort(ci().nodeset.begin(), ci().nodeset.end());
      
        dest[i]->rhs_lambda = lambda;
        dest[i]->rhs_mode = Split;
        dest[i]->rhs_nodes = {this, nullptr};
      
        dest[i]->syncInformation(lambda);
      }

      // Make sure we got all the nodes
      if(DEBUG_MODE)
        ci().solver.checkKeyIsGone(ci().nodeset.begin(), ci().nodeset.end(), ci().key);

      for(int i = 0; i < 2; ++i) {
        dest[i]->checkKeySynced();

        auto keys = ci().solver.getNeighborhoodKeySet
          (dest[i]->ci().nodeset.begin(), 
           dest[i]->ci().nodeset.end(), 
           dest[i]->ci().key);

        for(uint key : keys) {
          assert_notequal(key, ci().key);

          TVRegPathSegment *rps = rpsLookup(key);
          rps->ci().neighbors.insert(dest[i]);
          assert(rps != this);

          dest[i]->ci().neighbors.insert(rps);
        }
      }

      // Go through and remove this from the neighbors... 
      for(TVRegPathSegment* n_rps : ci().neighbors) {
        n_rps->ci().neighbors.erase(this);
      }

      assert(dest[0]->ci().neighbors.find(this) == dest[0]->ci().neighbors.end());
      assert(dest[1]->ci().neighbors.find(this) == dest[1]->ci().neighbors.end());

      // Now clean this one up
      lhs_lambda = lambda;
      lhs_mode = Split;
      lhs_nodes = dest;

      // For optimization 
      if(dest[0]->ci().nodeset.size() > dest[1]->ci().nodeset.size())
        swap(lhs_nodes[0], lhs_nodes[1]);

        if(PRINT_PATH_MESSAGES)
          cout << "Split  " << ci().key << " into " << lhs_nodes[0]->ci().key << " and " << lhs_nodes[1]->ci().key << "." << endl;
    } 

    ////////////////////////////////////////////////////////////////////////////////
    // JOINS

  public:

    void setupFromJoin(TVRegPathSegment *rps1, 
                       TVRegPathSegment *rps2,
                       dtype join_lambda) {

      // cout << "SETTIng up from JOIN!!!!!" << endl;

      const auto& ci1 = rps1->ci();
      const auto& ci2 = rps2->ci();

      if(DEBUG_MODE) {

        assert(ci2.neighbors.find(rps1) != ci2.neighbors.end());
        assert(ci1.neighbors.find(rps2) != ci1.neighbors.end());
        
        assert_equal(lambdaOfJoin(rps1, rps2), join_lambda);
        assert_equal(lambdaOfJoin(rps2, rps1), join_lambda);
      }

      // Pull in the nodes in the previous two segments
      assert(ci().nodeset.empty());
      ci().nodeset.resize(ci1.nodeset.size() + ci2.nodeset.size());
      merge(ci1.nodeset.begin(), ci1.nodeset.end(), 
            ci2.nodeset.begin(), ci2.nodeset.end(), 
            ci().nodeset.begin());

      // Set up the initial parts
      rhs_mode = Join;
      rhs_lambda = join_lambda;
      rhs_nodes = {rps1, rps2};
   
      // cout << "rps1: " 
      //      << rps1->get_r_AtLambda(join_lambda-1) << ", "
      //      << rps1->get_r_AtLambda(join_lambda) << ", "
      //      << rps1->get_r_AtLambda(join_lambda+1) << "]" << endl;

      // cout << "rps2: " 
      //      << rps2->get_r_AtLambda(join_lambda-1) << ", "
      //      << rps2->get_r_AtLambda(join_lambda) << ", "
      //      << rps2->get_r_AtLambda(join_lambda+1) << "]" << endl;

#ifndef NDEBUG
      {
        dtype r1m1 = rps1->get_r_AtLambda(join_lambda-1);
        dtype r1p1 = rps1->get_r_AtLambda(join_lambda+1);

        dtype r2m1 = rps2->get_r_AtLambda(join_lambda-1);
        dtype r2p1 = rps2->get_r_AtLambda(join_lambda+1);
      
        assert_leq(min(r1m1, r1p1), max(r2m1, r2p1));
        assert_leq(min(r2m1, r2p1), max(r1m1, r1p1));
      }
#endif

#ifndef NDEBUG      
      // In debug node, we need the old ones for checks when deactivating the joined nodes
      ci().neighbors.clear();
      ci().neighbors.insert(ci1.neighbors.begin(), ci1.neighbors.end());
      ci().neighbors.insert(ci2.neighbors.begin(), ci2.neighbors.end());
#else
      bool ci1_big = ci1.neighbors.size() > ci2.neighbors.size(); 

      auto& nb_big   = (ci1_big ? ci1.neighbors : ci2.neighbors);
      auto& nb_small = (!ci1_big ? ci1.neighbors : ci2.neighbors);

      ci().neighbors = std::move(nb_big);
      ci().neighbors.insert(nb_small.begin(), nb_small.end());
#endif
      
      ci().neighbors.erase(rps1);
      ci().neighbors.erase(rps2);

      // add in this key to all the other neighbors.
      for(TVRegPathSegment* nb : ci().neighbors) {
        auto& nbci = nb->ci();
        nbci.neighbors.insert(this);
        nbci.neighbors.erase(rps1);
        nbci.neighbors.erase(rps2);
      }

      if(DEBUG_MODE) {
        rps1->ci().neighbors.erase(rps2);
        rps2->ci().neighbors.erase(rps1);
      }

      // Clean up stuff in the joined segments 
      rps1->lhs_mode = Join;
      rps1->lhs_lambda = join_lambda;
      rps1->lhs_nodes = {this, 0};
      
      rps2->lhs_mode = Join;
      rps2->lhs_lambda = join_lambda;
      rps2->lhs_nodes = {this, 0};

      syncInformation(join_lambda);

      // Check a bunch of stuff
      
      assert_close(get_r_AtLambda(rhs_lambda), rps1->get_r_AtLambda(rps1->lhs_lambda), 1);
      assert_close(get_r_AtLambda(rhs_lambda), rps2->get_r_AtLambda(rps2->lhs_lambda), 1);

      checkKeySynced();

      if(PRINT_PATH_MESSAGES)
        cout << "Joined " << rps1->ci().key << " and " << rps2->ci().key << " into " << ci().key << "." << endl;
    }

    static inline dtype lambdaOfJoin(TVRegPathSegment* r1,
                                     TVRegPathSegment* r2) {
      return lambdaOfJoin(r1, r2, min(r1->rhs_lambda, r2->rhs_lambda));
    }

    static inline dtype lambdaOfJoin(TVRegPathSegment* r1,
                                     TVRegPathSegment* r2,
                                     dtype current_lambda) {
      r1->checkKeySynced();
      r2->checkKeySynced();

      // Make sure we are not joining one that's right here 
      if(r1->rhs_mode == Split && r2->rhs_mode == Split 
         && (r1->rhs_lambda == current_lambda && r2->rhs_lambda == current_lambda))
        return -1;

      dtype r10 = r1->adjusted_r_at_0;
      dtype r20 = r2->adjusted_r_at_0;
      dtype r11 = r1->adjusted_r_at_1;
      dtype r21 = r2->adjusted_r_at_1;
      
      dtype numer = r20 - r10;
      dtype denom = (r11 - r10) - (r21 - r20);

      dtype join_lambda = ( ( (numer < 0) != (denom < 0) )
                            ? -1
                            : Node::getScaleFromQuotient_T(numer, denom));

      // cout << "join_lambda = " << Node::scaleToValue(join_lambda) << endl;

      if(DEBUG_MODE) {
        assert_leq(current_lambda, r1->rhs_lambda);
        assert_leq(current_lambda, r2->rhs_lambda);
      
        // We've missed something if it's in this range

        if(join_lambda < min(r1->rhs_lambda, r2->rhs_lambda))
          assert_leq(join_lambda, current_lambda);
      }

      return (0 < join_lambda && join_lambda <= current_lambda) ? join_lambda : -1;
    }

    static inline dtype calculateJoins(TVRegPathSegment* r1,
                                       TVRegPathSegment* r2,
                                       dtype current_lambda) {

#ifndef NDEBUG
        const auto& ci1nb = r1->neighbors();
        const auto& ci2nb = r2->neighbors();

        assert(ci2nb.find(r1) != ci2nb.end());
        assert(ci1nb.find(r2) != ci1nb.end());
#endif      

      // Return the lambda at which these two segments join, or -1 if
      // they do not join at all.

      dtype join_lambda = lambdaOfJoin(r1, r2, current_lambda);

      if(join_lambda != -1) {
        r1->registerJoinPoint(join_lambda, r2);
        r2->registerJoinPoint(join_lambda, r1);
      }
       
      return join_lambda;
    }

    ////////////////////////////////////////////////////////////////////////////////
    // General node data

  public:
    size_t n_nodes;

    vector<node_ptr> nodeset;

    bool containsNode(node_ptr n) const {
      assert(!nodeset.empty());
      
      return binary_search(nodeset.begin(), nodeset.end(), n);
    }
      
    // This is set up 
    Array<TVRegPathSegment*, 2> rhs_nodes;
    Array<TVRegPathSegment*, 2> lhs_nodes;
  };

}; 


// #define NDEBUG
// #include "../common/debug.hpp"

#include "../common/debug_flymake_test.hpp"

#endif /* _TV_REGPATH_H_ */
