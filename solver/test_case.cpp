#include "solvers.hpp"

using namespace latticeQBP;

const static size_t size = 3;

int main(int argc, char **argv) {
  double f[9] = {0,1,2,3,4,5,6,7,8};
  vector<double> lambda(1);
  lambda[0] = 1;

  cout << "sizeof(long long) = " << sizeof(long long) << endl;

  calculate2dTV<Full2d_4, int64_t>(2,2,f,lambda);

  return 0;
}

