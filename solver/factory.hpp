#ifndef _FACTORY_H_
#define _FACTORY_H_

#include <exception>
#include <memory>
#include <string>

#include "interface.hpp"

namespace latticeQBP {

  typedef std::shared_ptr<LatticeFlowInterface> solver_ptr;

  using std::string;

  class BadSolverOrKernelName : public std::exception
  {
  public: 
    BadSolverOrKernelName(const string& _name) 
      : msg(string("Bad Kernel or Solver Name: ") + _name)
    {}

    virtual ~BadSolverOrKernelName() throw () {}

    virtual const char* what() const throw() 
    {
      return msg.c_str();
    }

  private:
    string msg;
  };

  // 1d interface
  solver_ptr getLatticeFlowInterface(
      const string& solver, const string& kernel, size_t nx);

  // 2d interface
  solver_ptr getLatticeFlowInterface(
				     const string& solver, const string& kernel, size_t nx, size_t ny);

  // 3d interface
  solver_ptr getLatticeFlowInterface(
				     const string& solver, const string& kernel, size_t nx, size_t ny, size_t nz);

  // 4d interface
  solver_ptr getLatticeFlowInterface(
				     const string& solver, const string& kernel, 
				     size_t nx1, size_t nx2, size_t nx3, size_t nx4);

  // Generic interface
  solver_ptr getLatticeFlowInterface(
				     const string& solver, const string& kernel, size_t *dimensions);

  bool isValidSolver(const string& _name);
  bool isValidKernel(const string& _name);

  int kernelDimension(const string& kernel);
  vector<string> validKernelNames();
  vector<string> validSolverNames();
}
#endif /* _FACTORY_H_ */
