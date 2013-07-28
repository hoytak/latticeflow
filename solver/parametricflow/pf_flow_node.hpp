#ifndef _PF_FLOW_NODE_H_
#define _PF_FLOW_NODE_H_

#include "../common.hpp"
#include "../network_flow/network_flow_node.hpp"

#include <algorithm>
#include <limits>

namespace latticeQBP {

  using namespace std;

// #ifdef NDEBUG
// #undef NDEBUG
// #endif
  
// #include "../common/debug.hpp"

  // The overall policy class 

  template <int _using_weights, int _using_scales, int _weights_binary>
  struct _PFNodePolicyBase {
    static constexpr bool on_by_reduction = true;
    static constexpr bool enable_key_partitioning = true;
    static constexpr bool adjustment_mode = true;
    static constexpr bool using_weights = _using_weights;
    static constexpr bool weights_binary = _weights_binary;
    static constexpr bool using_scales = _using_scales;
  };

  typedef _PFNodePolicyBase<false, false, false> PFUnweightedNodePolicy;
  typedef _PFNodePolicyBase<true, false, false>  PFWeightedNodePolicy;
  typedef _PFNodePolicyBase<true, false, true>  PFBinaryWeightedNodePolicy;

  typedef _PFNodePolicyBase<false, true, false> PFScaledUnweightedNodePolicy;
  typedef _PFNodePolicyBase<true, true, false>  PFScaledWeightedNodePolicy;
  typedef _PFNodePolicyBase<true, true, true>  PFScaledBinaryWeightedNodePolicy;

  ////////////////////////////////////////////////////////////////////////////////
  // Some stuff to deal with weights

  template <typename dtype, int using_weights, 
            int bits_precision, int weights_binary> struct __weightVariable {};

  template <typename dtype, int bits_precision, int weights_binary> 
  struct __weightVariable<dtype, 0, bits_precision, weights_binary> {

    static constexpr long mult_factor = 1;
    static inline void set(dtype v) { assert(false); }
    static inline dtype weight() { return 1; }
    static inline dtype mult(dtype x) { return x; }
    static inline dtype fvToNf(dtype x) { return x; }
    static inline dtype nfToFv(dtype x) { return x; }

  };

  template <typename dtype, int bits_precision, int weights_binary> 
  struct __weightVariable<dtype, 1, bits_precision, weights_binary> {
    __weightVariable() : _weight(0) {} 
    __weightVariable(double w) { set(w); } 
    
    static constexpr long mult_factor = weights_binary ? 1 : (long(1) << bits_precision);

    void set(double w) { 
      assert_leq(w, 1.0);
      assert_leq(0.0, w);
      _weight = dtype(round(w * mult_factor));
    }

    inline dtype weight() const { return _weight; }
    inline dtype mult(dtype x) const { return x * _weight; }
    inline dtype fvToNf(dtype x) const { return x * mult_factor; }
    inline dtype nfToFv(dtype x) const { return x / mult_factor; }

  private:
    dtype _weight;
  };
  
  ////////////////////////////////////////////////////////////////////////////////
  // The main class for this 

  template <class Kernel, typename dtype, typename NodePolicy> class PFFlowNode 
    : public NetworkFlowNode<Kernel, dtype, NodePolicy>
  {
    typedef NodePolicy Policy;
    typedef NetworkFlowNode<Kernel, dtype, Policy> Base;

  public:
    PFFlowNode() 
      :  base_fv(0)
      , current_fv(0)
      , current_fv_offset(0)
#ifndef NDEBUG
      , current_scale(-1)
#endif
    {}
    
    static_assert(!Policy::using_weights, "Weights not implemented yet.");

  public:
    // With the weights, we have that 
    static constexpr int n_weight_bits = (Policy::using_weights && !Policy::weights_binary
                                          ? 2*sizeof(dtype) : 0);
    
    static constexpr int n_scale_bits_room = 2*sizeof(dtype);
    static constexpr int n_bits_scale_precision = 7*sizeof(dtype) - n_scale_bits_room;

    static constexpr int log2_scale_max = n_bits_scale_precision + n_scale_bits_room;

    static constexpr int n_bits_function_room = (5*sizeof(dtype)) / 2;
    static constexpr int n_bits_function_precision = (7*sizeof(dtype) 
                                                      - n_weight_bits 
                                                      - n_bits_function_room);

    static constexpr int log2_function_max = n_bits_function_room + n_bits_function_precision;



  protected:
    __weightVariable<dtype, Policy::using_weights, 
                     n_weight_bits, Policy::weights_binary> _weight;
 
    dtype base_fv;
    dtype current_fv;
    dtype current_fv_offset;

#ifndef NDEBUG
    dtype current_scale;
#endif
 
  public:
    typedef typename CompType<dtype>::Type comp_type;
  
    static inline dtype toFVDType(double fv) { 
      double v1 = round(fv * double(dtype(1) << n_bits_function_precision));
      dtype v2 = dtype(v1);
      assert_equal(v1, double(v2));
      return v2;
    }

    static inline double toFValue(dtype fv) {
      return double(fv) / double(dtype(1) << n_bits_function_precision);
    }

    static inline dtype toScaleDType(double lm) {
      const int _n_bits_scale_precision = n_bits_scale_precision;

      assert(Policy::using_scales);
      assert_leq(0, lm);
      assert_leq(lm, double(dtype(1) << int(log2_scale_max)));

      return dtype(round(lm * double(dtype(1) << _n_bits_scale_precision)));
    }

    static inline double scaleToValue(dtype lm) {
      return double(lm) / double(dtype(1) << n_bits_scale_precision);
    }


    static inline dtype multFVScale(comp_type fv, dtype lm) {
      if(Policy::using_scales) {

        // assert_leq(abs(fv), (T(1) << log2_function_max));
        assert_leq(0, lm);
        assert_leq(lm, (dtype(1) << int(log2_scale_max)) );

        fv *= lm;
        fv += ((fv > 0) ? 1 : -1) * ( (dtype(1) << (n_bits_scale_precision - 1)));
        fv /= (dtype(1) << n_bits_scale_precision);
      }

      return toDType(fv);
    }

    template <typename T>
    static inline dtype getScaleFromQuotient_T(comp_type&& numer, const T& denom) {
      assert(Policy::using_scales);
      
      static constexpr dtype dtmax = numeric_limits<dtype>::max();

#ifndef NDEBUG
      comp_type old_numer = numer;
#endif
      numer *= (dtype(1) << n_bits_scale_precision);

      assert_equal(numer / (dtype(1) << n_bits_scale_precision), old_numer);

      numer /= denom;

#ifndef NDEBUG
      double v = static_cast<double>(old_numer) / static_cast<double>(denom);
      assert_close(v, static_cast<double>(numer) / (dtype(1) << n_bits_scale_precision), 1e-6);
#endif

      return likely(numer > dtmax) ? dtmax : toDType(numer);
    }

    template <typename T>
    static inline dtype getScaleFromQuotient(comp_type numer, const T& denom) {
      return getScaleFromQuotient_T(std::move(numer), denom);
    }

    template <typename T>
    static inline dtype castToScale(T scale) {

      assert_leq(0, scale);
      assert_leq(scale, (T(1) << int(log2_scale_max)));

      return toDType(scale);
    }

#ifndef NDEBUG
    dtype currentScale() const { return current_scale; }
#endif

    dtype fv() const {
      return base_fv;
    }

    dtype cfv() const {
      return current_fv;
    }

    dtype fv(dtype scale) const {
      assert(Policy::using_scales);
      return multFVScale(base_fv, scale);
    }

    dtype r() const {
      return Base::r() + current_fv_offset;
    }

    dtype qii() const {
      return base_fv;
    }

    dtype lm_qii() const {
      return current_fv;
    }

    dtype influence() const {
      return r() - lm_qii(); 
    }

    template <class Lattice> 
    void _debug_checkLevelsetMethodsNode(Lattice& lattice) {
#ifndef NDEBUG
      Base::_debugVerifyNodeConsistency(lattice);

      dtype cfo = current_fv_offset;
      dtype _cfv_predict = r();
      dtype _cfv = cfv();

      setOffset(lattice, 0);
      
      assert_equal(current_fv_offset, 0);
      assert_equal(r(), _cfv_predict);
      assert_equal(cfv(), _cfv);

      setOffset(lattice, 56323);
      
      assert_equal(current_fv_offset, 56323);
      assert_equal(r(), _cfv_predict);
      assert_equal(cfv(), _cfv);

      setOffset(lattice, cfo);

      assert_equal(current_fv_offset, cfo);
      assert_equal(r(), _cfv_predict);
      assert_equal(cfv(), _cfv);

      Base::_debugVerifyNodeConsistency(lattice);
#endif      
    }


    template <class Lattice> 
    inline void setOffsetAndScale(Lattice& lattice, dtype fv_offset, dtype scale)
    {
      assert(Policy::using_scales);

      dtype true_fv = fv(scale);

      dtype delta   = (true_fv - fv_offset) - (current_fv - current_fv_offset);

      _adjustValueInLattice(lattice, delta);

      current_fv        = true_fv;
      current_fv_offset = fv_offset;

#ifndef NDEBUG
      current_scale = scale;
#endif
    }

    template <class Lattice> 
    inline void setOffset(Lattice& lattice, dtype fv_offset) {
      dtype delta = (-fv_offset) - (- current_fv_offset);

      Base::_debugVerifyNodeConsistency(lattice);
      _adjustValueInLattice(lattice, delta);
      Base::_debugVerifyNodeConsistency(lattice);

      current_fv_offset = fv_offset;
    }

    template <class Lattice> 
    void _adjustValueInLattice(Lattice& lattice, dtype delta) {
      Base::adjustReduction(lattice, delta);
      //typename Base::template NodeFiller<Lattice>(lattice).addE1(this, 0, _weight.mult(delta)); 
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Convenience methods for working with the mean flow stuff.

    template <class Lattice, typename ForwardNodePtrIterator> 
    static inline bool setToMeanReductionPoint(Lattice& lattice,
                                               const ForwardNodePtrIterator& start, 
                                               const ForwardNodePtrIterator& end) {


      // Returns true if at the mean reduction; otherwise false.

      assert(start != end);

      comp_type r_sum = 0;
      double count = 0;
      dtype r_min = (*start)->r();
      dtype r_max = (*start)->r();

      // cout << "Reductions = ";

      // TO DO: implement this with weights
      for(ForwardNodePtrIterator it = start; it != end; ++it) {
        auto n = *it;
        n->_debug_checkLevelsetMethodsNode(lattice);
        dtype r = n->r();
        r_sum += r;
        r_min = min(r_min, r);
        r_max = max(r_max, r);
        count += 1;

        // cout << r << ",";
      }

      // cout << endl;

      // cout << "r_max = " << r_max << "; r_min = " << r_min << endl;

      if(r_max - r_min <= 1)
        return true;
      
      dtype r_new = dtype(ceil(double(r_sum) / count));

      // cout << "r_sum = " << r_sum << "; r_new = " << r_new;

      // To accomidate the case where flow still needs to be
      // redistributed, but it rounds to an end point.  
      if(r_new == r_max)
        --r_new;

      // cout << "; r_new_final = " << r_new << endl;

      // double base_level = 0; 

      // cout << "Levels: ";

      for(ForwardNodePtrIterator it = start; it != end; ++it) {
        auto n = *it;
        n->_debug_checkLevelsetMethodsNode(lattice);
        n->setOffset(lattice, r_new);
        n->_debug_checkLevelsetMethodsNode(lattice);
        // base_level += n->r();

        // cout << n->level() << ",";
      }
      // cout << endl;

      // cout << "base_level = " << (base_level / count) << endl;

      return false;
    }

    ////////////////////////////////////////////////////////////////////////////////
    // For initializing the lattice 

    void pullBaseFunctionValueFromLattice() {
      base_fv = Base::r();
      current_fv_offset = 0;
      current_fv = base_fv;

#ifndef NDEBUG
      current_scale = toScaleDType(1);
#endif

      assert_equal(r(), fv());
    }

    void setWeight(double w) {
      assert(Policy::using_weights);
      _weight.set(w);
    }

    // dtype level() const {
    //   return _weight.nfToFv(Base::level());
    // }
  };
};

// #define        NDEBUG
// #include "../common/debug.hpp"
  
#include "../common/debug.hpp"


#include "../common/debug_flymake_test.hpp"

#endif /* _TV_FLOW_NODE_H_ */
