#ifndef _TV_FLOW_NODE_H_
#define _TV_FLOW_NODE_H_

#include "../common.hpp"
#include "../network_flow/network_flow_node.hpp"

#include <algorithm>

namespace latticeQBP {

  using namespace std;

  template <class Kernel, typename dtype> class TVFlowNode 
    : public NetworkFlowNode<Kernel, dtype, NF_ENABLE_KEY_PARTIONING | NF_ADJUSTMENT_MODE> 
  {

    static constexpr int nf_parameters = NF_ENABLE_KEY_PARTIONING | NF_ADJUSTMENT_MODE;
    typedef NetworkFlowNode<Kernel, dtype, nf_parameters> Base;

  public:
    TVFlowNode() 
      : base_fv(0)
      , current_fv(0)
    {}
    
  protected:
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

      return (fv * lm) >> n_bits_lambda_precision;
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

    template <class Lattice>
    static void initTVLattice(const Lattice& lattice) {

      typedef typename Lattice::value_ptr node_ptr;

      typename Base::template NodeFiller<Lattice> filler(lattice);

      for(node_ptr n1 = lattice.begin(); n1 != lattice.end(); ++n1) {
        if(!lattice.withinBounds(n1))
          continue;

        filler.addE1(n1, 0, n1->fv());

        for(size_t ei = 0; ei < Lattice::Kernel_positive_size; ++ei) {
          node_ptr n2 = lattice.neighbor(n1, ei);

          if(lattice.withinBounds(n2)) {
            dtype v = dtype(round(lattice.geocutEdgeWeight(ei) * abs(n1->fv() - n2->fv())));
            filler.addE2(n1, n2, ei, 0, v, v, 0);
          }
        }
      }
    }

  };
};

#ifdef EMACS_FLYMAKE

#include "../kernels/kernels.hpp"

namespace latticeQBP {
  template class TVFlowNode<Star2d_4, long>;
};

#endif

#endif /* _TV_FLOW_NODE_H_ */
