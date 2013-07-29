#include "solvers.hpp"

using namespace latticeQBP;

const static size_t size = 3;
const static size_t nx = 10;
const static size_t ny = 10;

int main(int argc, char **argv) {
  double f[10] = {-1, 0,
                 0,  1, 
                 0.3,  1,  1,
                  -1,  0,  0};
  size_t nx = 2;
  size_t ny = 2;
  size_t nl = 10;

  vector<double> lambda(nl);

  for(size_t i = 0; i < lambda.size(); ++i)
    lambda[i] = 1000*double(i + 1) / lambda.size();

  auto R = calculate2dTV<Full2d_4, int64_t>(nx,ny,f,lambda);

  for(size_t xi = 0; xi < nx; ++xi) {
    for(size_t yi = 0; yi < ny; ++yi) {    
      cout << "R[" << xi << ',' << yi << "] = ";

      for(size_t i = 0; i < lambda.size(); ++i) {
        cout << R->at(i, xi,yi) << ", "; 
      }
      cout << endl;
    }
  }


  for(size_t xi = 0; xi < nx; ++xi) {
    for(size_t yi = 0; yi < ny; ++yi) {    
      cout << "CORRECT R[" << xi << ',' << yi << "] = ";
      for(size_t i = 0; i < lambda.size(); ++i) {
        auto R2 = calculate2dTV<Full2d_4, int64_t>(nx,ny,f,lambda[i]);
        cout << R2->at(xi, yi) << ",\t";
      }
      cout << endl;
    }
  }

  return 0;
}

