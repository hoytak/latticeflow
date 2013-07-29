
      // if(unlikely(calc_lambda < lambda_lb || calc_lambda > rhs_lambda)) {
      //   dtype lm_true = BisectionCheck(lambda_lb, rhs_lambda);

      //   cout << "Lambda Bizzaro on " << ci().key << ": calc_lambda = " << Node::scaleToValue(calc_lambda) 
      //        << "; lm_true = " << Node::scaleToValue(lm_true) << endl;
      //   cout << "          >>>> Dtype: calc_lambda = " << (calc_lambda) 
      //        << "; lm_true = " << (lm_true) << endl;

      //     cout 
      //        << "c = " << c << ".0\n" 
      //        << "g1 = " << g1 << ".0\n"
      //        << "g0 = " << g0 << ".0\n"
      //        << "q1 = " << q1 << ".0\n"
      //        << "q0 = " << q0 << ".07\n" 
      //        << "s1 = " << s1 << ".0\n"
      //        << "s0 = " << s0 << ".0\n" 
      //        << "lm = " << Node::scaleToValue(lm_true) << "\n" 
      //        << "dtlm = " << lm_true
      //        << "(s1*g0 - s0*g1) = " << toDType(s1*g0 - s0*g1) << ".0\n"
      //        << "c * (s0 + s1) = " << toDType(c * (s0 + s1)) << ".0\n"
      //        << "(s1*g0 - s0*g1) + c * (s0 + s1) = " << toDType((s1*g0 - s0*g1) + c * (s0 + s1)) << ".0\n"
      //        << "q1*s0 - q0*s1 = " << toDType(q1*s0 - q0*s1) << ".0\n"
      //        << endl;

      //   calc_lambda = lm_true;

      //   return DetailedSplitInfo({true, true, lm_true, cut_ptr});
      // } else
      {
      }

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

      assert_geq(calc_lambda, lambda_lb);

      // cout << "Calculated split lambda = " << calc_lambda << endl;

      // ci().solver.disableChecks();

      return DetailedSplitInfo({true, false, calc_lambda, cut_ptr});




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

#endif

      // Get the rest of the components to calculate the shape
      // comp_type qii_total   = p_info[0].qii_sum   + p_info[1].qii_sum;
      // comp_type gamma_total = p_info[0].gamma_sum + p_info[1].gamma_sum;

      // Now go through and see which one has the largest lambda 

      // typedef typename TV_PR_Class::PartitionInfo PartitionInfo;
      
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

      // auto predict = [&,s0,s1,g0,g1,c,q0,q1](int fg0, int fg1, int fc, int fq) {
      //   comp_type intercept = (fg0*s1*g0 - s0*fg1*g1) + fc*c*(s0 + s1);
      //   comp_type denom = q1*s0 - fq*q0*s1;
        
      //   dtype calc_lambda = Node::getScaleFromQuotient_T(std::move(intercept), denom);
        
      //   cout 
      //   << "Predicted lambda = " << calc_lambda << "(" << Node::scaleToValue(calc_lambda) << "), " 
      //   << fg0 << ":" << fg1 << ":" << fc << endl;
        
      //   return calc_lambda;
      // };

      // dtype lm_true = BisectionCheck(lambda_lb, rhs_lambda);
      // cout << "LM True = " << lm_true << endl;

      // cout 
      //   << "c = " << c << ".0\n" 
      //   << "g1 = " << g1 << ".0\n"
      //   << "g0 = " << g0 << ".0\n"
      //   << "q1 = " << q1 << ".0\n"
      //   << "q0 = " << q0 << ".07\n" 
      //   << "s1 = " << s1 << ".0\n"
      //   << "s0 = " << s0 << ".0\n" 
      //   << "lm = " << Node::scaleToValue(lm_true) << "\n" 
      //   << "dtlm = " << lm_true
      //   << endl;

      // predict(1,1,1, 1);
      // predict(1,1,1, -1);
      // predict(1,-1,1,1);
      // predict(1,-1,1,-1);
      // predict(1,-1,-1,1);
      // predict(1,-1,-1,-1);
      // predict(-1,1,1,1);
      // predict(-1,1,1,-1);
      // predict(-1,1,-1,1);
      // predict(-1,1,-1,-1);
      // predict(-1,-1,1,1);
      // predict(-1,-1,1,-1);
      // predict(-1,-1,-1,1);
      // predict(-1,-1,-1,-1);


    
    dtype findExactPointOfCutFeasibility(dtype pred_lambda, dtype lambda_lb, 
                                         typename TV_PR_Class::cutinfo_ptr cut) const {      

      cout << "Cut value = " << cut->cut_value << endl;

      auto score = [&](dtype lm) {
        // We're mainly concenrned about the point at which the positive excess 
        ci().solver.setRegionToLambda(ci().nodeset.begin(), ci().nodeset.end(), lm);
        
        dtype excess_0 = 0;
        dtype excess_1 = 0;
        for(node_ptr n : cut->partitions[1]->nodes) 
          excess_1 += n->current_excess();

        for(node_ptr n : cut->partitions[0]->nodes) 
          excess_0 += -n->current_excess();

        // cout << "excess at " << lm << " = " << excess_0 << "," << excess_1 << endl;

        return cut->cut_value - (excess_0 + excess_1);
      };

      cout << "pred_lambda = " << pred_lambda << ": " << score(pred_lambda) << endl;
      cout << "lambda_lb = " << lambda_lb << ": " << score(lambda_lb) << endl;
      cout << "rhs_lambda = " << rhs_lambda << ": " << score(rhs_lambda) << endl;

      // Now run the bisection algorithm to get the bounds to try and find the best version
      dtype lb, ub; 
      dtype score_lb, score_ub;
      
      assert_lt(score(lambda_lb), 0);
      assert_geq(score(rhs_lambda), 0);

      dtype plm_score = score(pred_lambda);
      bool searching_right;
      
      if(plm_score < 0) {
        lb = pred_lambda;
        score_lb = plm_score;
        searching_right = true;
      } else {
        ub = pred_lambda;
        score_ub = plm_score;
        searching_right = false;
      }

      // Now go back until it's negative 
      for(int shift_attempt = 8;;++shift_attempt) {
          
        dtype query_lm = pred_lambda + (searching_right ? 1 : -1) * (dtype(1) << shift_attempt);
        
        if(unlikely(!searching_right && query_lm <= lambda_lb)) {
          lb = lambda_lb; 
          score_lb = score(lb);
          assert_lt(score_lb, 0);
          break;
        }

        if(unlikely(searching_right && query_lm >= rhs_lambda)) {
          ub = rhs_lambda;
          score_ub = score(rhs_lambda);
          assert_geq(score_ub, 0);
          break;
        }

        dtype qs = score(query_lm);

        cout << "(" << lb << ":" << score(lb) << ") : (" << query_lm << "," << qs << ") : (" << ub << "," << score(ub) << ")" << endl;

        if(qs < 0) {
          lb = query_lm;
          score_lb = qs;
          if(!searching_right) 
            break;
        } else {
          ub = query_lm;
          score_ub = qs;
          if(searching_right)
            break;
        }
      }

      double weight_linear = 1;

      while(lb + 1 < ub) {
        cout << "lb = " << lb << "; ub = " << ub << endl;

        assert_lt(lb, ub);
        assert_equal(score_lb, score(lb));
        assert_equal(score_ub, score(ub));
        assert_lt(score_lb, 0);
        assert_geq(score_ub, 0);
        
        // The function is assumed to be close to linear.  Thus we
        // should start by attempting to nail that.   Error on the side 
        dtype mid_point = ((ub + lb) / 2);

        // double query_lm_dbl = (weight_linear * (lb + ((double(-score_lb) * (ub - lb)) / (score_ub - score_lb))) 
        //                        + (1 - weight_linear)*mid_point);

        // weight_linear *= 0.75;

        // assert_lt(query_lm_dbl, ub);
        // assert_gt(query_lm_dbl, lb); 

        // dtype query_lm = dtype( (query_lm_dbl < mid_point) ? ceil(query_lm_dbl) : floor(query_lm_dbl));

        dtype query_lm = mid_point;

        dtype qs = score(query_lm);

        if(qs < 0) {
          lb = query_lm;
          score_lb = qs;
        } else {
          ub = query_lm;
          score_ub = qs;
        }
      }

      return ub;
    }
