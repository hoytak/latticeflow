#include "direct_interface.hpp"
#include "../kernels/kernels.hpp"
#include "../common.hpp"
#include "../solvers.hpp"

#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace latticeQBP {
  using namespace std;

  ////////////////////////////////////////////////////////////////////////////////
  // This allows us to
%(instantiation_suppression_list)s

  ////////////////////////////////////////////////////////////////////////////////
  // The initialization functions

%(factory_creation_functions)s

   ////////////////////////////////////////////////////////////////////////////////
   // the lookup table

  typedef solver_ptr (*init_function)(size_t*);

  static const unordered_map<string, init_function> _generator_map = {
%(factory_generator_map)s
  };

  static const vector<string> _solver_list = {
%(factory_solver_list)s
  };

  static const vector<string> _kernel_list = {
%(factory_kernel_list)s
  };

  static const unordered_set<string> _kernel_set(_kernel_list.begin(), _kernel_list.end());
  static const unordered_set<string> _solver_set(_solver_list.begin(), _solver_list.end());

  static const unordered_map<string, int> _kernel_dimension = {
%(factory_kernel_dimension_list)s
  };

  ////////////////////////////////////////////////////////////////////////////////
  // 1d interface
  solver_ptr getLatticeFlowInterface(
      const string& solver, const string& kernel, size_t nx) {

    size_t buf[4] = {nx, 0, 0, 0};
    return getLatticeFlowInterface(solver, kernel, buf);
  }

  // 2d interface
  solver_ptr getLatticeFlowInterface(
      const string& solver, const string& kernel, size_t nx, size_t ny) {

    size_t buf[4] = {nx, ny, 0, 0};
    return getLatticeFlowInterface(solver, kernel, buf);
  }

  // 3d interface
  solver_ptr getLatticeFlowInterface(
      const string& solver, const string& kernel, size_t nx, size_t ny, size_t nz) {

    size_t buf[4] = {nx, ny, nz, 0};
    return getLatticeFlowInterface(solver, kernel, buf);
  }

  // 4d interface
  solver_ptr getLatticeFlowInterface(
     const string& solver, const string& kernel,
     size_t nx1, size_t nx2, size_t nx3, size_t nx4) {

    size_t buf[4] = {nx1, nx2, nx3, nx4};

    return getLatticeFlowInterface(solver, kernel, buf);
  }

  // Generic interface
  solver_ptr getLatticeFlowInterface(
     const string& solver, const string& kernel, size_t *dimensions) {

    string key = solver + ":" + kernel;

    auto it = _generator_map.find(key);

    if(it == _generator_map.end()) {
      throw BadSolverOrKernelName(solver + ", " + kernel);
    } else {
      return (it->second)(dimensions);
    }
  }

  bool isValidSolver(const std::string& _name) {
    return _solver_set.find(_name) != _solver_set.end();
  }

  bool isValidKernel(const std::string& _name) {
    return _kernel_set.find(_name) != _kernel_set.end();
  }

  int kernelDimension(const std::string& kernel) {

    auto it = _kernel_dimension.find(kernel);

    if(it == _kernel_dimension.end()) {
      throw BadSolverOrKernelName(kernel);
    } else {
      return (it->second);
    }
  }

  vector<string> validKernelNames() {
    return _kernel_list;
  }

  vector<string> validSolverNames() {
    return _solver_list;
  }

}
