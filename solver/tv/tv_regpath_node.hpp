#ifndef _TV_REGPATH_NODE_H_
#define _TV_REGPATH_NODE_H_

#include "../common.hpp"
#include "tv_flow_node.hpp"
#include "tv_push_relabel.hpp"

#include <set>
#include <map>
#include <vector>

namespace latticeQBP {
  
  using namespace std;

  template <typename dtype, typename TV_PR_Class> class TVRegPathSegment {
  public:
    typedef typename TV_PR_Class::KernelLattice Lattice;
    typedef typename TV_PR_Class::Node Node;    
    typedef typename TV_PR_Class::node_ptr node_ptr;    
    typedef typename CompType<dtype>::Type comp_type;
    

    typedef enum  {Unset, Initial, Join, Split} Mode;    

    TVRegPathSegment(uint _key, Lattice& _lattice, TV_PR_Class& _solver)
      : rhs_mode(Unset)
      , lhs_mode(Unset)
      , rhs_lambda(0)
      , lhs_lambda(DEBUG_MODE ? -1 : 0)
      , rhs_nodes(nullptr)
      , lhs_nodes(nullptr)
      , _construction_info(new ConstructionInfo(_key, _lattice, _solver))
    {}

    template <typename ForwardIterator>
    void setupAsInitial(const ForwardIterator& start, 
                        const ForwardIterator& end,
                        dtype lamba) {
      assert(nodeset.empty());

      nodeset.insert(start, end);
      sort(nodeset.start(), nodeset.end());

      syncKey(nodeset.start(), nodeset.end());
      
      rhs_lambda = lamba;
      rhs_mode = Initial;
      rhs_nodes = {nullptr, nullptr};

      updateLevelInformation(nodeset.start(), nodeset.end());
    }
    

    ////////////////////////////////////////////////////////////////////////////////
    // Internal data 
    
  private: 

    ////////////////////////////////////////////////////////////////////////////////
    // Stuff for tracking things 
    Mode rhs_mode, lhs_mode;
    dtype rhs_lambda, lhs_lambda;

    ////////////////////////////////////////////////////////////////////////////////
    // DIRECTIONAL INFORMATION for testing 

    // These are shifted by lambda_precision_bits to maintain some
    // additional numerical accuracy.
    dtype adjusted_r_at_1;
    dtype adjusted_r_at_0;

    template <typename T>
    static inline T adjust_r(T r) { return r * (T(1) << Node::n_bits_lambda_precision); }

    template <typename T>
    static inline T deadjust_r(T r) { return r / (T(1) << Node::n_bits_lambda_precision); }

    inline dtype get_r_AtLambda(dtype lambda) const {
      auto mult = [](dtype x, dtype lm){
        return Node::template multFVLambda<comp_type>(x, lm);};

      return dtype(deadjust_r(mult(adjusted_r_at_0, 1) 
                              + mult(adjusted_r_at_1 - adjusted_r_at_0, lambda)));
    }

    struct RegionInformation {
      
      // These are referenced from the base zero_reference.  
      dtype lm_qii_sum_d, qii_sum_d, r_sum_d;
      dtype zero_reference;
      dtype partition_size;

      comp_type qii_sum;

      // These are absolute values
      dtype gamma_sum;

    }; 

    template <typename ForwardIterator>
    inline RegionInformation getRegionInfo(const ForwardIterator& start, 
                                           const ForwardIterator& end, 
                                           dtype check_lambda
                                           ) {

      RegionInformation ri = {0,0,0,0,0,0};

      ri.zero_reference = (*start)->fv_predict();

      for(ForwardIterator it = start; it != end; ++it) {

        node_ptr n = *it;

        assert_equal(n->cfv(), n->fv(check_lambda));
        
        // Subtracting by fv_avg is for numerical stability
        ri.lm_qii_sum_d += n->cfv() - ri.zero_reference;
        ri.qii_sum_d += n->fv() - ri.zero_reference;

        ri.r_sum_d += n->cfv_predict() - ri.zero_reference;
        
        ++ri.partition_size;
      }
          
      ri.gamma_sum = (ri.r_sum_d - ri.lm_qii_sum_d);

      return ri;
    }

    template <typename ForwardIterator>
    void syncKey(const ForwardIterator& start, 
                 const ForwardIterator& end) {

      const auto& ci = *constructionInfo();

      n_nodes = 0;
      for(ForwardIterator it = start; it != end; ++it) {
        (*it)->setKey(ci.key);
        ++n_nodes;
      }
    }

    template <typename ForwardIterator>
    void updateLevelInformation(const ForwardIterator& start, 
                                const ForwardIterator& end, 
                                dtype check_lambda) {

      RegionInformation ri = getRegionInfo(start, end, check_lambda);

      // Get the information from these values to get the updated
      // directional information

      // The value of 
      adjusted_r_at_0 = dtype(adjuste_r(comp_type(ri.gamma_sum)) / ri.partition_size);
      
      comp_type r_sum_at_0 = comp_type(ri.zero_reference) * ri.partition_size + ri.qii_sum_d;

      adjusted_r_at_1 = dtype(adjuste_r(r_sum_at_0) / ri.partition_size);

      // Now, make sure that it's correct...

      if(DEBUG_MODE) {
        comp_type r_calc = get_r_AtLambda(check_lambda);

        dtype rhs_r = dtype( (comp_type(ri.zero_reference)*ri.partition_size + ri.r_sum_d) 
                             / ri.partition_size);

        assert_equal(r_calc, rhs_r);
      }
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
    ConstructionInfo* constructionInfo() const {
      assert(_construction_info != nullptr);
      return _construction_info;
    }

    set<TVRegPathSegment*>& neighbors() const {
      assert(_construction_info != nullptr);
      return _construction_info->neighbors;
    }

    void registerJoinPoint(dtype join_lambda, TVRegPathSegment* rps) const {
      constructionInfo()->join_points.insert(make_pair(join_lambda, rps));
    }

    dtype firstJoinPoint() const {
      assert(!constructionInfo()->join_points.empty());
      return constructionInfo()->join_points.rbegin()->first;
    }

    // Conditionally activate the nodes 
    void deactivate() {
      assert(lhs_lambda != -1);

      ConstructionInfo& ci = *constructionInfo();

      // Depending on the rhs mode, create a set with all the active
      // nodes in it. 
      if(rhs_mode == Split || rhs_mode == Initial) {
        nodeset.swap(ci.nodeset);
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
          assert(rhs_lambda > lhs_lambda);
          break;

        case Unset:
          assert(rhs_mode != Unset);
          break;
        }

        // Ensure it is not in any neighborhood maps
        for(TVRegPathSegment* rps : ci.neighbors) {
          assert(rps->constructionInfo()->neighbors.find(this) 
                 == rps->constructionInfo()->neighbors.end());
        }
      }

      delete constructionInfo();
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

      if(n_nodes == 1) {
        constructionInfo()->split_calculation_done_to_lambda = 0;
        constructionInfo()->lambda_of_split = -1;
        
        return SplitInfo({false, -1, 0});
      }
      
      auto& ci = (*constructionInfo());

      if(DEBUG_MODE && ci.split_calculation_done_to_lambda != -1)
        assert_leq(current_lambda, ci.split_calculation_done_to_lambda);

      // Now, go through and calculate a likely lower bound on the
      // split point.  There are situations where this ends up going
      // off, and will require recalculations, but those are few and
      // far between.
      
      dtype lambda_calc_lb = 0;

      for(const auto& p : ci.join_points) {
        if(p.first == p.second->firstJoinPoint()) {
          lambda_calc_lb = max(p.first, lambda_calc_lb);
        }
      }
      
      DetailedSplitInfo dsi, last_dsi;
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

        if(DEBUG_MODE) {
          auto dsi2 = _calculateSingleSplit(dsi.lambda_of_split_capacity - 1);
          assert_equal(dsi2.lambda_of_split_capacity, dsi.lambda_of_split_capacity);
        }

        lambda_calc = dsi.lambda_of_split_capacity;
      }

      // Store that information in the computation structure
      if(cut_exists) {
        assert(last_dsi.split_occurs);
        
        ci.split_information = last_dsi.cut;

      } else {
        ci.split_calculation_done_to_lambda = lambda_calc_lb;
        ci.lambda_of_split = -1;
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

    DetailedSplitInfo _calculateSingleSplit(dtype lambda_lb) const {

      const ConstructionInfo& ci = *constructionInfo();

      // Calculate the split point...
      ci.solver.setRegionToLambda(ci.nodeset.begin(), 
                                  ci.nodeset.end(), lambda_lb);
      
      // Now see if it can be solved at that lambda
      ci.solver.runSection(ci.nodeset.begin(), ci.nodeset.end(), ci.key);
      auto cut_ptr = ci.solver.runPartitionedSection(ci.nodeset.begin(), ci.nodeset.end(), ci.key);
      const auto& piv = cut_ptr->partitions;

      assert_equal(piv.size(), 2);

      if(!cut_ptr->any_cut) 
        return DetailedSplitInfo({false, 0, nullptr});
      
      Array<RegionInformation, 2> p_info = {
        getRegionInfo(piv[0].begin(), piv[0].end(), lambda_lb),
        getRegionInfo(piv[1].begin(), piv[1].end(), lambda_lb),
      };

      // Get the rest of the components to calculate the shape
      comp_type qii_total   = p_info[0].qii_sum   + p_info[1].qii_sum;
      comp_type gamma_total = p_info[0].gamma_sum + p_info[1].gamma_sum;

      // Now go through and see which one has the largest lambda 
      dtype max_lambda_so_far = 0;

      typedef typename TV_PR_Class::PartitionInfo PartitionInfo;

      auto calcLambda = [&, qii_total, gamma_total](const RegionInformation& ri, const PartitionInfo& pi) {
        size_t R_size = ri.partition_size;

        comp_type lambda_coeff = R_size * ri.qii_sum   - ri.size * qii_total;
        comp_type lambda_intcp = R_size * comp_type(ri.gamma_sum) - ri.size * gamma_total;
        comp_type cut = R_size * comp_type(ri.pt->cut_value);

        // is_on being true means that this partition was on the
        // high end, so at the lambda = 0 end, it's got too much
        // flow if this is the blocking cut section.  This means
        // that the incoming flow must decrease with increasing
        // lambda, and that the original intercept term must be
        // positive.  Thus we are looking for the point where it
        // hits the cut.

        assert(! ((ri.pt->is_on  && ( lambda_coeff >= 0 || cut >= lambda_intcp)  )
                  || (!ri.pt->is_on && ( lambda_coeff <= 0 || cut >= -lambda_intcp) ) ));
            
        dtype calc_lambda = Node::getLambdaFromQuotient(abs(lambda_intcp) - cut, abs(lambda_coeff));

        assert_geq(calc_lambda, 0);
        assert_leq(calc_lambda, rhs_lambda);

        return calc_lambda;
      };

      dtype calc_lambda = calcLambda(p_info[0], piv[0]);
      
      assert_equal(calc_lambda, calcLambda(p_info[1], piv[1]));

      return DetailedSplitInfo({true, calc_lambda, cut_ptr});
    }


  public:
    template <typename ForwardIterator, typename RPSLookupFunction>
    void setupFromSplit(const ForwardIterator& start, 
                        const ForwardIterator& end,
                        TVRegPathSegment* parent,
                        dtype lambda, const RPSLookupFunction& rpsLookup) {

      auto& ci = *constructionInfo();
      
      assert(ci.nodeset.empty());

      ci.nodeset.insert(start, end);
      sort(ci.nodeset.start(), ci.nodeset.end());
      
      syncKey(ci.nodeset.start(), ci.nodeset.end());
      
      rhs_lambda = lambda;
      rhs_mode = Split;
      rhs_nodes = {parent, nullptr};
      
      updateLevelInformation(ci.nodeset.start(), ci.nodeset.end());

      // Build the neighborhood map.
      ci.solver.constructNeighborhoodSet(ci.nodeset.start(), ci.nodeset.end(), ci.key, 
                                         [&](node_ptr n) {
                                           TVRegPathSegment *nrps = rpsLookup(n->key());
                                           ci.neighbors.insert(nrps);
                                         });
    }

    template 
    void applySplit(TVRegPathSegment *dest1, TVRegPathSegment *dest2, dtype lambda) {
      
      const auto& ci = *constructionInfo();

      assert_equal(lambda, ci.lambda_of_split);
      
      // First, ensure that the cut is applied to the section so the
      // appropriate edges are saturated.
      ci.solver.applyPartioningCut(ci.cut_ptr, ci.key);
      
      // Split up the nodes
      const auto& piv = ci.cut_ptr->partitions;

      dest1->setupAsSplit(piv[0].start(), piv[0].end(), lambda);
      dest2->setupAsSplit(piv[1].start(), piv[1].end(), lambda);
      
      // Now clean this one up
      lhs_lambda = lambda;
      lhs_mode = Split;
      lhs_nodes = {dest1, dest2};
    } 

    ////////////////////////////////////////////////////////////////////////////////
    // JOINS

  public:

    void setupFromJoin(dtype join_lambda, 
                       TVRegPathSegment *rps1, 
                       TVRegPathSegment *rps2) {

      assert(rps2->neighbor_keys.find(rps1->key) != rps2->neighbor_keys.end());
      assert(rps1->neighbor_keys.find(rps2->key) != rps1->neighbor_keys.end());

      assert_equal(lambdaOfJoin(rps1, rps2, false), join_lambda);
      assert_equal(lambdaOfJoin(rps2, rps1, false), join_lambda);

      // Pull in the nodes in the previous two segments
      assert(nodeset.empty());
      nodeset.resize(rps1->nodeset.size() + rps2->nodeset.size());
      merge(rps1->nodeset.begin(), rps1->nodeset.end(), 
            rps2->nodeset.begin(), rps2->nodeset.end(), 
            nodeset.begin());

      syncKey(nodeset.start(), nodeset.end());

      // Set up the initial parts
      rhs_mode = Join;
      rhs_lambda = join_lambda;
      rhs_nodes = {rps1, rps2};

      assert_equal(rps1->get_r_AtLambda(join_lambda), 
                   rps2->get_r_AtLambda(join_lambda));

      ConstructionInfo& ci = *constructionInfo();
      ConstructionInfo& rps1_ci = *(rps1->constructionInfo());
      ConstructionInfo& rps2_ci = *(rps2->constructionInfo());

      // Get the neighbors together.
      ci.neighbors = rps1_ci.neighbors + rps2_ci.neighbors;
      
      // ci.neighbors.insert(rps2->constructionInfo()->neighbors.begin(), 
      //                  rps2->constructionInfo()->neighbors.end());

      ci.neighbors.erase(rps1);
      ci.neighbors.erase(rps2);

      // add in this key to all the other neighbors.
      for(TVRegPathSegment* nb : ci.neighbors) {\
        auto& nbci = *(nb->constructionInfo());
        nbci.neighbors.insert(this);
        nbci.neighbors.erase(rps1);
        nbci.neighbors.erase(rps2);
      }

      // Clean up stuff in the joined segments 
      rps1->lhs_mode = Join;
      rps1->lhs_lambda = join_lambda;
      rps1->lhs_nodes = {this, 0};
      
      rps2->lhs_mode = Join;
      rps2->lhs_lambda = join_lambda;
      rps2->lhs_nodes = {this, 0};

      updateLevelInformation(nodeset.start(), nodeset.end());

      // Check a bunch of stuff
      assert_equal(get_r_AtLambda(rhs_lambda), rps1->get_r_AtLambda(rps1->lhs_lambda));
      assert_equal(get_r_AtLambda(rhs_lambda), rps2->get_r_AtLambda(rps1->lhs_lambda));

    }

    static inline dtype calculateJoins(TVRegPathSegment* r1,
                                       TVRegPathSegment* r2,
                                       dtype current_lambda) {

      if(DEBUG_MODE) {
        assert(r1->neighbors.find(r2) != r1->neighbors.end());
        assert(r2->neighbors.find(r1) != r2->neighbors.end());
      }

      // Return the lambda at which these two segments join, or -1 if
      // they do not join at all.

      // Do the fast check first
      if( (r1->adjusted_r_at_0 < r2->adjusted_r_at_0) 
          == (r1->adjusted_r_at_r1 < r2->adjusted_r_at_r1) )
        return -1;
      
      comp_type r10 = r1->adjusted_r_at_0;
      comp_type r20 = r2->adjusted_r_at_0;
      comp_type r11 = r1->adjusted_r_at_1;
      comp_type r21 = r2->adjusted_r_at_1;

      comp_type adj_join_lambda = getLambdaFromQuotient(r11 - r10, (r20 - r10) - (r21 - r11));
      dtype join_lambda = Node::castToLambda(deadjust_r(adj_join_lambda));
      
      assert(join_lambda != 0);
      assert_leq(join_lambda, r1->rhs_lambda);
      assert_leq(join_lambda, r2->rhs_lambda);

      if(0 < join_lambda && join_lambda < current_lambda) {
        r1->registerJoinPoint(join_lambda, r2);
        r2->registerJoinPoint(join_lambda, r1);
        return join_lambda;
      } else {
        return -1;
      }
    }



    ////////////////////////////////////////////////////////////////////////////////
    // General node data
  private:


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
      
    size_t n_nodes;

    vector<node_ptr> nodeset;
      
    // This is set up 
    Array<TVRegPathSegment*, 2> rhs_nodes;
    Array<TVRegPathSegment*, 2> lhs_nodes;
  };

}; 

#ifdef EMACS_FLYMAKE

#include "../kernels/kernels.hpp"

namespace latticeQBP {
  template class TVRegPathSegment<long, _TV_PRFlow_Test>;
};
#endif

#endif /* _TV_REGPATH_H_ */
