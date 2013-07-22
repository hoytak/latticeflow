#include "solvers.hpp"

using namespace latticeQBP;

const static size_t size = 3;

int main(int argc, char **argv) {
  double f[9] = {0,1,-1,
                 0,1,1,
                 -1,0,0};

  vector<double> lambda(1);
  lambda[0] = 1;

  calculate2dTV<Full2d_4, int64_t>(3,3,f,lambda);

  return 0;
}

