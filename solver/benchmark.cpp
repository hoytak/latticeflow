#ifdef ENABLE_OPENMP
#undef ENABLE_OPENMP
#endif

#include "graphcuts/bk_energy.hpp"

#include "energy.hpp"
#include "rng.hpp"
#include "indexing.hpp"
#include "utilities.hpp"
#include "timetracker.hpp"
#include "array_wrapper.hpp"
#include "kernels/kernels.hpp"

typedef long dtype;

using namespace latticeQBP;

//#define _LatticeEnergy NetFlowLatticeEnergyReductions
#define _LatticeEnergy LatticeEnergyMinimizer

#include <boost/lexical_cast.hpp>
#include <string>
#include <iostream>
#include <algorithm>
#include <boost/algorithm/string.hpp> 

using namespace std;
using namespace boost;

struct ImageOptions {
  ImageOptions () 
    : noise_lb(0)
    , noise_ub(100000)
    , n_points(30)
    , points_lb(-50000)
    , points_ub(50000)
  {}
  int noise_lb, noise_ub;
  int n_points;
  int points_lb, points_ub;
  int maxDiff() const {
    int v = (1 + noise_ub - noise_lb) + (1 + points_ub - points_lb);
    assert_geq(v, 0);
    return v;
  }
};

// Not very fast, but okay for now
template <int n_dimensions>
pair<size_t, int> argclosest(const vector< Array<size_t, n_dimensions> >& pts, 
			     const Array<size_t, n_dimensions>& x) {
  size_t idx = 0;
  int min_value = -1;

  for(size_t i = 0; i != pts.size(); ++i) {
    int d = 0;
    for(size_t j = 0; j != pts[i].size(); ++j) {
      d += (pts[i][j] - x[j]) * (pts[i][j] - x[j]);
    }

    if( d > min_value) {
      idx = i;
      min_value = d;
    }
  }

  return make_pair(idx, min_value);
}

template <int nd> Array<size_t, nd> coords(size_t idx, size_t edge_size) {

  Array<size_t, nd> c;

  for(int j = nd; j != 0; --j) {
    c[j-1] = (idx % edge_size);
    idx /= edge_size;
  }

  return c;
}

template <int nd>
static void construct_random_bitmap(vector<dtype>& X, size_t edge_size, 
				    const ImageOptions& opt, unsigned long seed) {
  
  FastInt rng(seed); 

  // First build the correct array
  Array<size_t, nd> dims;
  fill(dims.begin(), dims.end(), edge_size);
  Array<size_t, nd> dim_product(dim_vector_factors<nd>(dims));
  
  // Construct the random points
  vector< Array<size_t, nd> > points(opt.n_points);
  vector< int > point_values(opt.n_points);
  
  for(size_t i = 0; i < points.size(); ++i) {
    for(size_t j = 0; j < nd; ++j)
      points[i][j] = rng(edge_size);
    
    point_values[i] = rng(opt.points_lb, opt.points_ub);
  }

  // Setup the random matrix
  X.resize(prod<size_t, nd>(dims));
  Array<size_t, nd> x;

  for(IndexIterator<nd> idxit(dims); !idxit.done(); ++idxit) {
    size_t pidx = argclosest<nd>(points, idxit.coords()).first;
    X[idxit.index()] = point_values[pidx] + rng(opt.noise_lb, opt.noise_ub);
  }

  cout << "Image constructed; " << X.size() << " pixels." << endl;
}

template <int nd> class GraphCutEnergyWrapper {
public:
  typedef typename Indexer<nd>::index_vect index_vect;

  GraphCutEnergyWrapper(const index_vect& _dimensions) 
    : dimensions(_dimensions)
    , indexer(_dimensions)
    , graphcutter(indexer.size(), nd*indexer.size())
    , graph_nodes(indexer.size())
  {
    for(size_t i = 0; i < indexer.size(); ++i) {
      graph_nodes[i] = graphcutter.add_variable();
    }
  }

  void addE1(const index_vect& coord, long E0, long E1) {
    graphcutter.add_term1(graph_nodes[indexer[coord]], E0, E1);
  }

  void addE2(const index_vect& c1, 
	     const index_vect& c2,
	     int E00, int E01, int E10, int E11) {
    
    graphcutter.add_term2(graph_nodes[indexer[c1]], 
			  graph_nodes[indexer[c2]],
			  E00, E01, E10, E11);
  }

  void run() {
    graphcutter.minimize();
  }

  bool on(const index_vect& node_index) {
    return graphcutter.get_var(graph_nodes[indexer[node_index]]);
  }
  
private:
  index_vect dimensions; 
  Indexer<nd> indexer;
  typedef Energy<dtype, dtype, dtype> GraphCutter;
  GraphCutter graphcutter;
  vector<typename GraphCutter::Var> graph_nodes;
};

template <int nd, typename LatticeEnergy> 
void make(LatticeEnergy *le,
	  GraphCutEnergyWrapper<nd> *gc,
	  char mode, const ImageOptions& opt, 
	  const vector<dtype>& X, size_t edge_size) {

  typedef typename Indexer<nd>::index_vect index_vect;

  index_vect dimensions(edge_size); 
  Indexer<nd> indexer(dimensions);

  if(mode == 's' || mode == 'S') {
    
    for(auto it = le->getPairwiseFillingIterator(); !it.done(); ++it) {

      const index_vect& idx_1 = it.latticeCoordOf1();
      const index_vect& idx_2 = it.latticeCoordOf2();

	double d = (double((X[indexer[idx_1]] - X[indexer[idx_2]])) / opt.maxDiff() );
	int sameness = int(ceil((65536.0 / dist2(idx_1, idx_2)) 
				* exp(-d*d) ));

	gc->addE2(idx_1, idx_2, -sameness, 0, 0, -sameness);
	it.addPairwisePotential(-sameness, 0, 0, -sameness);
    }
   
    // Now fix some of the points to zero or one; need to fix a swath
    // so it has an effect.
    FastInt rng(size_t( (500 + X[0] + X[1]) )); 
    
    for(size_t i = 0; i < 10; ++i) {
      int s = i % 2;
      int energy_value = 65536 * 100;
      
      size_t index = rng(X.size());
      int dim = rng(nd);
      int dir = rng(2) ? -1 : 1;

      for(size_t j = 0; j < edge_size / 5; ++j) {
	size_t set_index = indexer.neighborIndex(index, dim, dir*j);
	if(set_index >= X.size()) {
	  break;
	} else {
	  le->addUnaryPotential(indexer[set_index], -s*energy_value, -(1-s)*energy_value);
	  gc->addE1(indexer[set_index], -s*energy_value, -(1-s)*energy_value);
	}
      }
    }

    cout << "Problem constructed." << endl; 

  } else {
    
    cout << "Mode " << mode << " not recognized." << endl;
    exit(1);
  };
}

template <int nd, typename Kernel> void run_benchmark(char mode, size_t edge_size, size_t random_seed) {

  typedef typename Indexer<nd>::index_vect index_vect;
  
  vector<dtype> X;

  ImageOptions opt;

  construct_random_bitmap<nd>(X, edge_size, opt, 10000 + random_seed);

  index_vect dimensions(edge_size); 

  // Now run it...
  TimeTracker tt;

  tt.start();
  _LatticeEnergy<nd, Kernel, dtype> le_qbp(dimensions);
  GraphCutEnergyWrapper<nd> le_gc(dimensions);

  make<nd>(&le_qbp, &le_gc, mode, opt, X, edge_size);

  tt.stop();

  cout << "Time taken in setup = " << tt.asString() << "." << endl;


  cout << "Starting LatticeEnergy version." << endl;
  tt.reset();
  tt.start();
  le_qbp.run();
  tt.stop();
  cout << "Time taken in QBP optimization = " << tt.asString() << "." << endl;

  cout << "Starting GraphCuts version." << endl;
  tt.reset();
  tt.start();
  le_gc.run();
  tt.stop();
  cout << "Time taken in GC optimization = " << tt.asString() << "." << endl;

  size_t n_pos = 0;
  size_t n_neg = 0;
  size_t fpos_mismatch_count = 0;
  size_t fneg_mismatch_count = 0;
  
  for(IndexIterator<nd> idxit(dimensions); !idxit.done(); ++idxit) {
    bool q_on  = le_qbp.on(idxit.coords());
    bool gc_on = le_gc.on(idxit.coords());

    if(!q_on && gc_on)
      ++fneg_mismatch_count;

    if(q_on && !gc_on)
      ++fpos_mismatch_count;

    ++(q_on ? n_pos : n_neg);
  }
  
  if( (fneg_mismatch_count + fpos_mismatch_count) != 0) {
    cout << "WARNING: mismatch count between QBP and GC is " 
	 << (fneg_mismatch_count + fpos_mismatch_count) 
	 << " (false negatives = " << fneg_mismatch_count
	 << ", false positives = " << fpos_mismatch_count 
	 << ", total negatives = " << n_neg
	 << ", total positives = " << n_pos
	 << ")"
	 << "!" << endl;
  } else {
    cout << "Solutions exactly match." << endl;
  }
}


int main(int argn, char** argv)
{ 
  if(argn != 4 && argn != 5) {
    cout << "\nUsage: " << argv[0] << " <mode> <n dimensions> <edge size> [n_times]"  
	 << endl;
    return 1;
  }

  // Okay
  char mode = argv[1][0];
  size_t n_dim = lexical_cast<size_t>(argv[2]);
  size_t edge_size = lexical_cast<size_t>(argv[3]);
  size_t n_times = (argn >= 5) ? lexical_cast<size_t>(argv[4]) : 1;

  for(size_t nt = 0; nt < n_times; ++nt) {
    switch(n_dim) {

    case 2:
      run_benchmark<2, Full2d_4>(mode, edge_size, nt);
      break;
    case 21:
      run_benchmark<2, Full2d_12>(mode, edge_size, nt);
      break;
    case 22:
      run_benchmark<2, Full2d_48>(mode, edge_size, nt);
      break;
    case 3:
      run_benchmark<3, Full3d_6>(mode, edge_size, nt);
       break;
    case 31:
      run_benchmark<3, Full3d_32>(mode, edge_size, nt);
       break;

    default:
      cerr << "\nBenchmarks not supported for diminsion " << n_dim << "." << endl;
      return 1;
    }
  }
  
  return 0;
}


