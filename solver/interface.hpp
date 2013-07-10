#ifndef _INTEFACE_H_
#define _INTEFACE_H_

#include "common.hpp"

#include <exception>
#include <memory>
#include <string>
#include <vector>


namespace latticeQBP {
  using std::vector;
  using std::string;

  class LatticeFlowInterface {
  
  public:

    ////////////////////////////////////////////////////////////////////////////////
    // Stuff for interfacing easily to external interfaces without
    // dealing with the pain of the templating features.  Note that this
    // class has to be created through the factory functions.

    // The dimension of the lattice.
    size_t nDimensions() const { return _n_dimensions; }

    // Returns the per-dimension sizing of the lattice.
    const vector<size_t>& shape() const { return _shape; }

    // Total number of nodes in the lattice. Equal to the product of the
    // sizes returned by shape().
    size_t nNodes() const { return _n_nodes; }

    // Returns the number of edges to set per node, as determined by the
    // kernel.  This is equal to half the number of edges connected to a
    // typical node -- thus each edge is considered once in the 
    size_t nFillingEdges() const { return _n_filling_edges; }
  
    // Returns a vector of size nFillingEdges() * nDimensions(), in
    // which entry [i * nDimensions() + j] gives the offset of edge i in
    // dimension j.  
    const vector<int>& edgeDeltas() const {
      return _edge_deltas;
    }
  
    //////////////////////////////////////////////////////////////////////////////////
    // Filling the lattice values.  
    // 
    // The unary_energy_array must be of length nNodes()*2, where
    // unary_energy_array[2*node_index + 0] gives the potential for node
    // node_index being off and unary_energy_array[2*node_index + 1]
    // gives the potential for node_index being on.
    //
    // The binary_energy_array must be of length
    // ``nNodes()*nFillingEdges()*4``, where entry
    // ``binary_energy_array[node_index*(4*nFillingEdges()) +
    // edge_index*4 + potential_index]`` gives the energy for edge
    // edge_index from node node_index and potential potential_index.
    // The different potentials are set in order E00, E01, E10, E11. 

    virtual void addEnergyPotentials(long *unary_energy_array, long *binary_energy_array) = 0;

    ////////////////////////////////////////////////////////////////////////////////
    // Conversion routines to make the indexing easier


    // Converts an array of length nDimensions() to the absolute
    // node_index used in filling and retrieving values
    inline size_t nodeIndex(size_t *node_coords) const {
      size_t mult = 1;
      size_t node_index = 0;

      for(size_t dim = _n_dimensions; (dim--) != 0;) {
        node_index += mult * node_coords[dim];
        mult *= _shape[dim];
      }
      return node_index;
    }

    // Returns the index in unary_energy_array corresponding to the off
    // potential of node node_index,
    inline size_t unary_E0_index(size_t node_index) const {
      return 2*node_index + 0;
    }

    // Returns the index in unary_energy_array corresponding to the on
    // potential of node node_index,
    inline size_t unary_E1_index(size_t node_index) const {
      return 2*node_index + 1;
    }

    inline size_t binary_E00_index(size_t node_index, size_t edge_index) const {
      return node_index * (4*_n_filling_edges) + 4*edge_index + 0;
    }

    inline size_t binary_E01_index(size_t node_index, size_t edge_index) const {
      return node_index * (4*_n_filling_edges) + 4*edge_index + 1;
    }

    inline size_t binary_E10_index(size_t node_index, size_t edge_index) const {
      return node_index * (4*_n_filling_edges) + 4*edge_index + 2;
    }

    inline size_t binary_E11_index(size_t node_index, size_t edge_index) const {
      return node_index * (4*_n_filling_edges) + 4*edge_index + 3;
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Control functions 

    virtual void run() = 0;

    virtual double runTime() = 0;

    ////////////////////////////////////////////////////////////////////////////////
    // Functions to retrieve the results of the optimization.  For
    // simple cut problems, these will be just 0 or 1; for the
    // parametric flow problems or regularization methods.

    virtual vector<int> getCut() const = 0;

  protected:

    // To force creation only through the factory function
    LatticeFlowInterface(){}

    void init(size_t n_dimensions, 
              size_t n_filling_edges,
              vector<size_t> shape,
              vector<int> edge_deltas
              )
    {
      _n_dimensions = n_dimensions;
      _n_filling_edges = n_filling_edges;
      _n_nodes = prod(shape);
      _shape = shape;
      assert_equal(edge_deltas.size(), _n_filling_edges * _n_dimensions);
      _edge_deltas = edge_deltas;
    }

  private:
    size_t _n_dimensions, _n_filling_edges, _n_nodes;
    vector<size_t> _shape;
    vector<int> _edge_deltas;
  };

  typedef std::shared_ptr<LatticeFlowInterface> solver_ptr;

  ////////////////////////////////////////////////////////////////////////////////
  // Generic creation functions 

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
  solver_ptr getLatticeFlowInterface(const string& solver, const string& kernel, size_t nx);

  // 2d interface
  solver_ptr getLatticeFlowInterface(const string& solver, const string& kernel, size_t nx, size_t ny);

  // 3d interface
  solver_ptr getLatticeFlowInterface(const string& solver, const string& kernel, size_t nx, size_t ny, size_t nz);

  // 4d interface
  solver_ptr getLatticeFlowInterface(const string& solver, const string& kernel, 
				     size_t nx1, size_t nx2, size_t nx3, size_t nx4);

  // Generic interface
  solver_ptr getLatticeFlowInterface(const string& solver, const string& kernel, size_t *dimensions);

  bool isValidSolver(const string& _name);
  bool isValidKernel(const string& _name);

  int kernelDimension(const string& kernel);
  vector<string> validKernelNames();
  vector<string> validSolverNames();
};


#endif /* _INTEFACE_H_ */
