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
    typedef typename TV_PR_Class::Lattice Lattice;
    typedef typename TV_PR_Class::Node Node;    
    typedef typename TV_PR_Class::node_ptr node_ptr;    
    typedef typename CompType<dtype>::comp_type comp_type;
    

    typedef enum  {Unset, Initial, Join, Split} Mode;    

    TVRegPathSegment(uint _key, Lattice& _lattice, TV_PR_Class& _solver)
      : key(_key)
      , lattice(_lattice)
      , solver(_solver)
      , rhs_mode(Unset)
      , lhs_mode(Unset)
      , split_point_calculated_to_lambda(-1)
    {}

    template <typename ForwardIterator, typename RPSKeyLookupFunction> 
    void setupAsInitial(const ForwardIterator& start, 
                        const ForwardIterator& end,
                        dtype solved_lamba,
                        const RPSKeyLookupFunction& key_lookup
                        ) {

      for(ForwardIterator it = start; it != end; ++it)
        (*it)->setKey(key);
      
      assert(nodeset.empty());
      nodeset.insert(start, end);
      
      rhs_lambda = solved_lamba;
      rhs_mode = Initial;

      updateLevelInformation(start, end);

      // Add up the 

    }

    ////////////////////////////////////////////////////////////////////////////////
    // Internal data 
    
  private: 

    const uint key;
    Lattice& lattice;
    TV_PR_Class& solver;

    ////////////////////////////////////////////////////////////////////////////////
    // Stuff for tracking things 
    Mode rhs_mode;

    dtype rhs_lambda;
    dtype rhs_r;

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

      rhs_r = dtype( (comp_type(ri.zero_reference)*ri.partition_size + ri.r_sum_d) 
                     / ri.partition_size);

      if(DEBUG_MODE) {
        comp_type r_calc = Node::multFVLambda(adjusted_r_at_1 - adjusted_r_at_0, 
                                              comp_type(check_lambda));

        assert_equal(rhs_r, r_calc);
      }
    }

    ////////////////////////////////////////////////////////////////////////////////
    // SPLITS

  public:

    struct SplitInfo {
      TVRegPathSegment* splitting_node;
      
      vector<node_ptr> split_low, split_high;
      
      vector<pair<node_ptr, uint> > cut_edges;
    };
    
    typedef shared_ptr<SplitInfo> splitinfo_ptr;

    struct SplitPointInfo {
      dtype split_lambda;
      splitinfo_ptr info;
    };

    splitinfo_ptr getSplitPoint(dtype lambda_lb) {

      if(unlikely(split_point_calculated_to_lambda != -1 
                  && lambda_lb >= split_point_calculated_to_lambda) ) {
        return spi;
      }

      // Calculate the split point...
      
      
      // Store that information in here...
      
      return spi;
    }

  private:
    dtype split_point_calculated_to_lambda;

    SplitPointInfo spi;

    ////////////////////////////////////////////////////////////////////////////////
    // JOINS

  public:

    static void join(dtype join_lambda, 
                     TVRegPathSegment* dest,
                     TVRegPathSegment *rps1, 
                     TVRegPathSegment *rps2) {

      assert(rps2->neighbor_keys.find(rps1->key) != rps2->neighbor_keys.end());
      assert(rps1->neighbor_keys.find(rps2->key) != rps1->neighbor_keys.end());

      assert_equal(lambdaOfJoin(rps1, rps2, false), join_lambda);
      assert_equal(lambdaOfJoin(rps2, rps1, false), join_lambda);

      dest->rhs_mode = Join;
      dest->rhs_lambda = join_lambda;
      dest->rhs_nodes = {rps1, rps2};
      dest->rhs_r = rps1->get_r_AtLambda(join_lambda);

      assert_equal(rps1->get_r_AtLambda(join_lambda), 
                   rps2->get_r_AtLambda(join_lambda));

      // Get the neighbors together.
      dest->neighbor_keys = rps1->neighbor_keys;
      dest->neighbor_keys.insert(rps2->neighbor_keys.begin(), 
                                 rps2->neighbor_keys.end());

      dest->neighbor_keys.erase(rps1->key);
      dest->neighbor_keys.erase(rps2->key);

      // add in this key to all the other neighbors.
      for(const auto& neighbor_key_ptr : dest->neighbor_keys) {
        uint key = neighbor_key_ptr.first;
        TVRegPathSegment* nb = neighbor_key_ptr.second;
        
        nb->neighbor_keys.insert(make_pair(dest->key, dest));
        nb->neighbor_keys.erase(rps1->key);
        nb->neighbor_keys.erase(rps2->key);
      }

      // Clean up stuff in the joined segments 
      rps1->lhs_mode = Join;
      rps1->lhs_lambda = join_lambda;
      rps1->lhs_nodes = {dest, 0};
      
      rps2->lhs_mode = Join;
      rps2->lhs_lambda = join_lambda;
      rps2->lhs_nodes = {dest, 0};
    }

    static inline dtype lambdaOfJoin(TVRegPathSegment* r1,
                                     TVRegPathSegment* r2,
                                     bool neighbor_check) {

      // Return the lambda at which these two segments join, or -1 if
      // they do not join at all.

      // Do the fast check first
      if( (r1->adjusted_r_at_0 < r2->adjusted_r_at_0) 
          == (r1->adjusted_r_at_r1 < r2->adjusted_r_at_r1) )
        return -1;

      if(neighbor_check) {
        if(r1->neighbor_keys.find(r2->key) == r1->neighbor_keys.end())
          return -1;
      }
      
      comp_type r10 = r1->adjusted_r_at_0;
      comp_type r20 = r2->adjusted_r_at_0;
      comp_type r11 = r1->adjusted_r_at_1;
      comp_type r21 = r2->adjusted_r_at_1;

      comp_type adj_join_lambda = getLambdaFromQuotient(r11 - r10, (r20 - r10) - (r21 - r11));
      dtype join_lambda = Node::castToLambda(deadjust_r(adj_join_lambda));
      
      assert(join_lambda != 0);
      assert_leq(join_lambda, r1->rhs_lambda);
      assert_leq(join_lambda, r2->rhs_lambda);

      return dtype(join_lambda);
    }

  private:

    ////////////////////////////////////////////////////////////////////////////////
    // General data


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
      
    vector<node_ptr> nodeset;
      
    map<uint, TVRegPathSegment*> neighbor_keys;

    // This is set up 
    Array<TVRegPathSegment*, 2> rhs_nodes;
    Array<TVRegPathSegment*, 2> lhs_nodes;
  };

}; 

#endif /* _TV_REGPATH_H_ */
