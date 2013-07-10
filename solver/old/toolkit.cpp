#include "array_wrapper.hpp"
#include "indexing.hpp"
#include "energy.hpp"
#include "toolkit.hpp"
#include <cmath>

using namespace latticeQBP;

/*

template <int nd> void pinnedLevelSets(double *L, double *X, const Array<size_t, nd>& x_nd, 
				       const Array<double, nd>& weights) {
  
  typedef typename Indexer<nd>::index_vect index_vect;

  NetFlowLatticeEnergyReductions<nd> e(x_nd);

  cout << "Setting up problem with dimensions " << x_nd << endl;

  Indexer<nd> indexer(x_nd);

  e.fillE2([&, X, indexer, weights](const index_vect& idx_1, const index_vect& idx_2) {
      
      Array<double, nd> correlations(0);
      
      double diff = (X[indexer[idx_1]] - X[indexer[idx_2]]);

      // double tie_value = exp( - diff*diff);
      
      double tie_value = 1 / (0.1 + abs(diff));// * (weights * abs(idx_2 - idx_1)).sum();

      return Array<double, 4>({-tie_value, 0, 0, -tie_value});
    });


  for(auto it = IndexIterator<nd-1>(slice<1,nd>(x_nd)); !it.done(); ++it) {
    e.addE1(cat(0, it.coords()), -1, 0);
    e.addE1(cat(x_nd[0]-1, it.coords()), 0, -1);
  }

  // 
  e.run();

  for(auto it = IndexIterator<nd>(x_nd); !it.done(); ++it)
    L[it.index()] = e.level(it.coords());
}

void pinnedLevelSets(double *L, double *X, 
		     size_t nz, size_t ny, size_t nx, 
		     double wx, double wy, double wz) {

  pinnedLevelSets(L, X, Array<size_t, 3>({nz, ny, nx}), Array<double, 3>({wz, wy, wx}));
}
			   
////////////////////////////////////////////////////////////////////////////////

template <int nd, typename Kernel, int n_covariates> 
void getImageLevelSets(double *L, double *X, double *predictions, double *labels,
		       const Array<size_t, nd>& dimensions, double e2_balance, double label_value) {

  typedef typename Indexer<nd>::index_vect index_vect;

  auto to_dtype  = [](double x) { return long(round(65536.0 * x) ); };
  auto to_double = [](long x)   { return double(x) / 65536.0; };

  NetFlowLatticeEnergyReductions<nd, Kernel, long> e(dimensions);

  // get the average variance of the kernel covariates 
  Array<vector<double>, n_covariates> values;
  vector<double> weights;

  for(auto& v : values)
    v.reserve(Kernel::positive_size * prod(dimensions));

  weights.reserve(Kernel::positive_size * prod(dimensions));

  // Now, go through and fill these

  Indexer<nd+1> indexer(cat(dimensions, n_covariates));
  Indexer<nd> pr_indexer(dimensions);

  e.dryRunE2([&, X, indexer, &values, &weights, to_dtype]
	     (const index_vect& idx_1, const index_vect& idx_2) {

	       double d = 1.0 / dist(idx_1, idx_2);
	       weights.push_back(d);
      
	       for(int i = 0; i < n_covariates; ++i) {
		 double diff = X[indexer[cat(idx_1, i)]] - X[indexer[cat(idx_2, i)]];
		 values[i].push_back(diff);
	       }
	     }
	     );
  
  // Calculate the mean and variance
  Array<double, n_covariates> means, variances;
  
  for(int i = 0; i < n_covariates; ++i) {
    means[i] = inner_product(values[i].begin(), values[i].end(), 
			     weights.begin(), double(0)) / values[i].size();
    
    variances[i] = 0;

    for(size_t j = 0; j < values[i].size(); ++j) {
      double v = values[i][j] - means[i];
      variances[i] += v*v*weights[j];
    }
      
    variances[i] /= values[i].size();
  }
  
  // Get the norm
  double norm = 1;//accumulate(weights.begin(), weights.end(), double(0)) / prod(dimensions);

  // Now, go through and run with it
  size_t idx = 0;
  e.fill(
	 [&, predictions, labels, label_value, pr_indexer, e2_balance, to_dtype](const index_vect& idx_1) {
	   size_t idx = pr_indexer[idx_1];
	   double v = labels[idx] * label_value + (1 -e2_balance) * predictions[idx];

	   return Ar(long(0), to_dtype(-v)); 
	 },

	 [&, &idx, values, variances, weights, norm]
	 (const index_vect& idx_1, const index_vect& idx_2) {
	   
	   double v = 0;
	   
	   for(int i = 0; i < n_covariates; ++i) {
	     v += values[i][idx] * values[i][idx] / (2 * variances[i] * norm);
	   }

	   ++idx;

	   long lv = to_dtype(e2_balance * weights[idx] * exp(-v));

	   return Ar(long(0), lv, lv, long(0));
	 }
	 );
 
  // Okay, now that we know it, run it!
  cout << "Model built, now running it." << endl;

  e.run();

  for(auto it = IndexIterator<nd>(dimensions); !it.done(); ++it)
    L[it.index()] = to_double(e.level(it.coords()));

}

void getRGBImageLevelSets(double *L, double *X, double *predictions, 
			  double *labels, 
			  size_t nx, size_t ny, double e2_balance, double label_value, 
			  int kernel) {

    // Labels is 0 for no label, 1 for label -> 1, and -1 for label 0. 
  switch(kernel) {
  case 0:
    return getImageLevelSets<2, Neighbors<2>, 3>(L, X, predictions, labels, Ar(nx, ny), e2_balance, label_value);
  case 1:
    return getImageLevelSets<2, DiagNeighbors<2>, 3>(L, X, predictions, labels, 
						     Ar(nx, ny), e2_balance, label_value);
  case 2:
    return getImageLevelSets<2, Star<2,2>, 3>(L, X, predictions, labels, Ar(nx, ny), e2_balance, label_value);
  case 3:
    return getImageLevelSets<2, Star<2,2>, 3>(L, X, predictions, labels, Ar(nx, ny), e2_balance, label_value);
  case 4:
    return getImageLevelSets<2, Star<2,2>, 3>(L, X, predictions, labels, Ar(nx, ny), e2_balance, label_value);
  case 5:
    return getImageLevelSets<2, Star<2,2>, 3>(L, X, predictions, labels, Ar(nx, ny), e2_balance, label_value);
  case 6:
    return getImageLevelSets<2, Star<2,2>, 3>(L, X, predictions, labels, Ar(nx, ny), e2_balance, label_value);

  default:;
  };
}

*/
