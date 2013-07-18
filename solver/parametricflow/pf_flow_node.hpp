#ifndef _PF_FLOW_NODE_H_
#define _PF_FLOW_NODE_H_

#include "../common.hpp"
#include "../network_flow/network_flow_node.hpp"

#include <algorithm>

namespace latticeQBP {

  using namespace std;

  static constexpr int _parflow_nf_parameters = NF_ENABLE_KEY_PARTIONING | NF_ADJUSTMENT_MODE;

  template <typename dtype, int using_weights> struct __weightVariable {};

  template <typename dtype> 
  struct __weightVariable<dtype, 0> {
    void setWeight(dtype v) { assert(false); }
    dtype getWeight() { return 1; }
  };

  template <typename dtype> 
  struct __weightVariable<dtype, 1> {
    __weightVariable() : _weight() {}

    void setWeight(double v) { _weight = v;}
    dtype getWeight() { return _weight; }
  private:
    dtype _weight;
  };

  
  ////////////////////////////////////////////////////////////////////////////////
  // The main class for this 

  template <class Kernel, typename dtype, int weighted> class PFFlowNode 
    : public NetworkFlowNode<Kernel, dtype, _parflow_nf_parameters>
  {
    static constexpr int nf_parameters = _parflow_nf_parameters;

    typedef NetworkFlowNode<Kernel, dtype, nf_parameters> Base;

  public:
    PFFlowNode() 
      : _weight(0)
      ,  base_fv(0)
      , current_fv(0)
      , current_fv_offset(0)
    {}
    
  protected:
    dtype _weight; 
    dtype base_fv;
    dtype current_fv;
    dtype current_fv_offset;
 
  public:
    typedef typename CompType<dtype>::Type comp_type;

    static int constexpr n_bits_function_room = 8;
    static int constexpr n_bits_lambda_max = sizeof(dtype)*1;

    static int constexpr n_bits_function_precision = (sizeof(dtype)*5 - n_bits_function_room);
    static int constexpr n_bits_lambda_precision = (sizeof(dtype)*3 - n_bits_lambda_max);
  
    static inline dtype toFVDType(double fv) {
      assert_leq(abs(fv), 1);
      return dtype(round(fv * double(dtype(1) << n_bits_function_precision)));
    }

    static inline double toFValue(dtype fv) {
      return double(fv) / double(dtype(1) << n_bits_function_precision);
    }

    static inline dtype toLmDType(double lm) {
      assert_geq(0, lm);
      assert_leq(lm, double(dtype(1) << n_bits_lambda_max));

      return dtype(round(lm * double(dtype(1) << n_bits_lambda_precision)));
    }

    template <typename T>
    static inline T multFVLambda(T fv, T lm) {
      assert_leq(abs(fv), (T(1) << (n_bits_function_precision + n_bits_function_room) ) );
      assert_geq(lm, 0);
      assert_leq(lm, (T(1) << (n_bits_lambda_precision + n_bits_lambda_max) ));

      // Rounds up to nearest value; This ensures that 0 <==> 0.0
      return (fv * lm + ( (T(1) << (n_bits_lambda_precision)) - 1) ) >> n_bits_lambda_precision;
    }

    template <typename T>
    static inline T getLambdaFromQuotient(T numer, T denom) {
      numer *= (T(1) << n_bits_lambda_precision);
      numer /= denom;
      
      return numer;
    }

    template <typename T>
    static inline dtype castToLambda(T lambda) {
      assert_geq(lambda, 0);
      assert_leq(lambda, (T(1) << (n_bits_lambda_precision + n_bits_lambda_max) ));
      return dtype(lambda);
    }

    dtype fv() const {
      return base_fv;
    }

    dtype cfv() const {
      return current_fv;
    }

    dtype fv(dtype lambda) {
      return multFVLambda(base_fv, lambda);
    }

    dtype cfv_predict() const {
      return Base::level() + current_fv_offset;
    }

    template <class Lattice> 
    void setFunctionValue(const Lattice& lattice, dtype fv_offset, dtype lambda) {
      dtype true_fv = fv(lambda);
      dtype delta   = (true_fv - fv_offset) - (current_fv - current_fv_offset);

      typename Base::template NodeFiller<Lattice>(lattice).addE1(this, 0, delta);

      current_fv        = true_fv;
      current_fv_offset = fv_offset;
    }

    template <class Lattice> 
    void setOffset(const Lattice& lattice, dtype fv_offset) {
      dtype delta = (-fv_offset) - (- current_fv_offset);

      typename Base::template NodeFiller<Lattice>(lattice).addE1(this, 0, delta); 

      current_fv_offset = fv_offset;
    }

    ////////////////////////////////////////////////////////////////////////////////
    // For initializing the lattice 

    void setBaseFunctionValue(double x) {
      base_fv = toFVDType(x);
    }
  };
};

#ifdef EMACS_FLYMAKE

#include "../kernels/kernels.hpp"

namespace latticeQBP {
  template class PFFlowNode<Star2d_4, long, true>;
  template class PFFlowNode<Star2d_4, long, false>;
};

#endif

#endif /* _TV_FLOW_NODE_H_ */
