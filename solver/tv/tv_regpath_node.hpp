#ifndef _TV_REGPATH_NODE_H_
#define _TV_REGPATH_NODE_H_

#include "../common.hpp"
#include "tv_push_relabel.hpp"

#include <set>
#include <map>
#include <vector>

namespace latticeQBP {

#ifdef NDEBUG
#undef NDEBUG
#endif

#include "../common/debug.hpp"


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
      checkKeySynced();

      assert(lhs_lambda != -1);

      // Depending on the rhs mode, create a set with all the active
      // nodes in it. 
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
#endif
      }

      assert(_construction_info != nullptr);
      delete _construction_info;
      _construction_info = nullptr;
    }


    ////////////////////////////////////////////////////////////////////////////////
    // SPLITS

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

      for(const auto& p : ci().join_points) {
        if(p.first == p.second->firstJoinPoint()) {
          lambda_calc_lb = max(p.first, lambda_calc_lb);
        }
      }
      
      DetailedSplitInfo dsi{false,0}, last_dsi;
      dtype lambda_calc = lambda_calc_lb;
      bool cut_exists = false;

      while(true) {
        last_dsi = dsi;
        dsi = _calculateSingleSplit(lambda_calc);

        if(!dsi.split_occurs)
          break;
        else
          cut_exists = true;

        assert_gt(dsi.lambda_of_split_capacity, lambda_calc);
        assert_leq(dsi.lambda_of_split_capacity, current_lambda);

        lambda_calc = dsi.lambda_of_split_capacity + 1;
      }

      // Store that information in the computation structure
      if(cut_exists) {

        assert(last_dsi.split_occurs);

        // if(DEBUG_MODE) {
        //   auto dsi2 = _calculateSingleSplit(last_dsi.lambda_of_split_capacity - 1);
        //   assert_equal(dsi2.lambda_of_split_capacity, last_dsi.lambda_of_split_capacity);
          
        //   if(dsi.like_trains_passing_in_the_dead_of_a_cold_moonless_night)
        // } 
        
        ci().split_information = last_dsi.cut;
        ci().lambda_of_split = lambda_calc;

        return SplitInfo({true, lambda_calc, lambda_calc_lb});

      } else {
        ci().split_calculation_done_to_lambda = lambda_calc_lb;
        ci().lambda_of_split = -1;
        return SplitInfo({false, -1, lambda_calc_lb});
      }
    }

  private:


    struct DetailedSplitInfo{ 
      bool split_occurs;

      // This gives the exact lambda of this cut; if it is one less
      // than this, then the cut will form.  Above this, there may be
      // a new cut, but at least this one won't form.
      dtype lambda_of_split_capacity;

      typename TV_PR_Class::cutinfo_ptr cut;
    };

    dtype BisectionCheck(dtype lambda_lb, dtype lambda_ub) const {

      while(lambda_lb + 1 < lambda_ub) {

        dtype try_lambda = (lambda_ub + lambda_lb) / 2;

        // Calculate the split point...
        ci().solver.setRegionToLambda(ci().nodeset.begin(), 
                                      ci().nodeset.end(), try_lambda, false);
      
        auto cut_ptr = 
          ci().solver.runPartitionedSection(ci().nodeset.begin(), ci().nodeset.end(), ci().key);


        if(!cut_ptr->any_cut)
          lambda_ub = try_lambda;
        else
          lambda_lb = try_lambda;
        
        // cout << "[" 
        //      << lambda_lb << ", " 
        //      << try_lambda << ", " 
        //      << lambda_ub << "]" << (cut_ptr->any_cut ? "Cut! " : "No Cut!") << endl;

      }

      cout << "A cut occurs at " << lambda_lb << " but not at " << lambda_ub << endl;
      return lambda_lb;
    }

    DetailedSplitInfo _calculateSingleSplit(const dtype lambda_lb) const {

      // ci().solver.enableChecks();

      cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
      cout << "Calculating split at lambda = " << lambda_lb << endl;

      // Calculate the split point...
      ci().solver.setRegionToLambda(ci().nodeset.begin(), 
                                    ci().nodeset.end(), lambda_lb, false);
      
      // Now see if it can be solved at that lambda
      auto cut_ptr = ci().solver.runPartitionedSection(ci().nodeset.begin(), ci().nodeset.end(), ci().key);
      const auto& piv = cut_ptr->partitions;

      if(!cut_ptr->any_cut) {
        cout << "     > No split. " << endl;
        ci().solver.disableChecks();
        return DetailedSplitInfo({false, 0, nullptr});
      } else {
        assert(!piv[0]->nodes.empty());
        assert(!piv[1]->nodes.empty());
      }

      Array<RegionInformation, 2> p_info = {
        getRegionInfo(piv[0]->nodes.begin(), piv[0]->nodes.end(), lambda_lb),
        getRegionInfo(piv[1]->nodes.begin(), piv[1]->nodes.end(), lambda_lb),
      };

      comp_type c = cut_ptr->cut_value;

      comp_type g0 = p_info[0].gamma_sum - c;
      comp_type g1 = p_info[1].gamma_sum + c;

      comp_type q0 = p_info[0].qii_sum;
      comp_type q1 = p_info[1].qii_sum;

      comp_type s0 = cut_ptr->partitions[0]->nodes.size();
      comp_type s1 = cut_ptr->partitions[1]->nodes.size();

      // comp_type g0 = cut_ptr->gamma_sum[0];
      // comp_type g1 = cut_ptr->gamma_sum[1];

      // comp_type q0 = cut_ptr->qii_sum[0];
      // comp_type q1 = cut_ptr->qii_sum[1];

      // comp_type s0 = cut_ptr->partitions[0]->nodes.size();
      // comp_type s1 = cut_ptr->partitions[1]->nodes.size();

      // comp_type c = cut_ptr->cut_value;


#ifndef NDEBUG

      RegionInformation ri_ref = getRegionInfo(ci().nodeset.begin(), ci().nodeset.end(), lambda_lb);

      assert_equal(ri_ref.gamma_sum, p_info[0].gamma_sum + p_info[1].gamma_sum);
      assert_equal(ri_ref.qii_sum,   p_info[0].qii_sum   + p_info[1].qii_sum);

      // dtype lm_true = BisectionCheck(lambda_lb, rhs_lambda);

      // cout << "cut_ptr->cut_value = " << cut_ptr->cut_value
      //      << "; p_info[1].gamma_sum = " << p_info[1].gamma_sum 
      //      << "; p_info[0].gamma_sum = " << p_info[0].gamma_sum
      //      << "; p_info[1].qii_sum = " << p_info[1].qii_sum 
      //      << "; p_info[0].qii_sum = " << p_info[0].qii_sum
      //      << endl;

      // cout 
      //   << "c = " << c << ".0\n" 
      //   << "g1 = " << g1 << ".0\n"
      //   << "g0 = " << g0 << ".0\n"
      //   << "q1 = " << q1 << ".0\n"
      //   << "q0 = " << q0 << ".07\n" 
      //   << "s1 = " << s1 << ".0\n"
      //   << "s0 = " << s0 << ".0\n" 
      //   << "lm = " << Node::scaleToValue(lm_true)
      //   << "dtlm = " << lm_true
      //   << endl;

#endif

      // Get the rest of the components to calculate the shape
      // comp_type qii_total   = p_info[0].qii_sum   + p_info[1].qii_sum;
      // comp_type gamma_total = p_info[0].gamma_sum + p_info[1].gamma_sum;

      // Now go through and see which one has the largest lambda 

      typedef typename TV_PR_Class::PartitionInfo PartitionInfo;
      
      // auto predict = [g0,g1,s0,s1,q0,q1,c](int f0, int f1, int f2) {
      //   comp_type g1p = g1*f1;
      //   comp_type g0p = g0*f0;
      //   comp_type cp = c*f2;

      //   comp_type intercept = (s1*g0p - s0*g1p) - cp*(s0 + s1);
      //   comp_type denom = q1*s0 - q0*s1;

      //   dtype calc_lambda = Node::getScaleFromQuotient_T(std::move(intercept), denom);

      //   cout << "Predicted lambda = " << calc_lambda << "(" << Node::scaleToValue(calc_lambda) << "), " 
      //   << f0 << ":" << f1 << ":" << f2 << endl;

      //   return calc_lambda;
      // }; 

      comp_type intercept = (s1*g0 - s0*g1) + c*(s0 + s1);
      comp_type denom = q1*s0 - q0*s1;

      dtype calc_lambda = Node::getScaleFromQuotient_T(std::move(intercept), denom);

      // dtype calc_lambda = predict(1,1,1);
      // predict(1,1,-1);
      // predict(1,-1,1);
      // predict(1,-1,-1);
      // predict(-1,1,1);
      // predict(-1,1,-1);
      // predict(-1,-1,1);
      // predict(-1,-1,-1);

      // is_on being true means that this partition was on the
      // high end, so at the lambda = 0 end, it's got too much
      // flow if this is the blocking cut section.  This means
      // that the incoming flow must decrease with increasing
      // lambda, and that the original intercept term must be
      // positive.  Thus we are looking for the point where it
      // hits the cut.

      // assert(! ((ri.pt->is_on  && ( lambda_coeff >= 0 || cut >= lambda_intcp)  )
      //           || (!ri.pt->is_on && ( lambda_coeff <= 0 || cut >= -lambda_intcp) ) ));

      assert_leq(calc_lambda, rhs_lambda);

      // if(cut_ptr->cut_value == 0)
      //   return DetailedSplitInfo({true, true, calc_lambda, cut_ptr});

      // assert_geq(calc_lambda, lambda_lb);

      cout << "Calculated split lambda = " << calc_lambda << endl;

      // ci().solver.disableChecks();

      return DetailedSplitInfo({true, calc_lambda, cut_ptr});
    }


  public:

    template <typename RPSLookupFunction>
    void applySplit(Array<TVRegPathSegment*, 2> dest,
                    dtype lambda, 
                    const RPSLookupFunction& rpsLookup) {
      
      assert_equal(lambda, ci().lambda_of_split);
      
      // First, ensure that the cut is applied to the section so the
      // appropriate edges are saturated.
      ci().solver.applyPartioningCut(ci().split_information, ci().key);
      
      // Split up the nodes
      const auto& piv = ci().split_information->partitions;

      for(int i = 0; i < 2; ++i) {

        assert(dest[i]->ci().nodeset.empty());

        dest[i]->ci().nodeset = move(ci().split_information->partitions[i]->nodes);
        // [i]->sort(ci().nodeset.begin(), ci().nodeset.end());
      
        dest[i]->rhs_lambda = lambda;
        dest[i]->rhs_mode = Split;
        dest[i]->rhs_nodes = {this, nullptr};
      
        dest[i]->syncInformation(lambda);


        dest[i]->checkKeySynced();
      }


      for(int i = 0; i < 2; ++i) {
        auto& nb = dest[i]->ci().neighbors;

        auto keys = ci().solver.getNeighborhoodKeySet
          (ci().nodeset.begin(), ci().nodeset.end(), ci().key);

        for(uint key : keys) {
          nb.insert(rpsLookup(key));
        }
      }

      // Now clean this one up
      lhs_lambda = lambda;
      lhs_mode = Split;
      lhs_nodes = dest;

      // For optimization 
      if(dest[0]->ci().nodeset.size() > dest[1]->ci().nodeset.size())
        swap(lhs_nodes[0], lhs_nodes[1]);

      // Go through and remove this from the neighbors... Man, I am
      // drunk on R. Jelinek fernet and a fine 7 yr. panama rum. 1:3
      // they make a fine manhattan but damn I need some bitterman's
      // burlesque bitters with it.  

      for(TVRegPathSegment* n_rps : ci().neighbors) {
        n_rps->ci().neighbors.erase(this);
      }
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

#define NDEBUG
#include "../common/debug.hpp"

#include "../common/debug_flymake_test.hpp"

#endif /* _TV_REGPATH_H_ */
