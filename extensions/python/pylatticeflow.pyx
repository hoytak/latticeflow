from libcpp.string cimport string
from libcpp.vector cimport vector
import numpy as np
import numpy.linalg as la
from numpy import dot, inner
from numpy cimport ndarray as ar
from cython cimport view
from cPickle import dumps
from hashlib import md5
cimport cython

cdef extern from "math.h":
    double fabs "abs"(double)
    double fsqrt "sqrt"(double)

cdef extern from "interface.hpp"  namespace "latticeQBP":

    cdef cppclass LatticeFlowInterface:

        size_t nDimensions()
        size_t nNodes()
        vector[size_t] shape()
        size_t nFillingEdges()
        vector[int] edgeDeltas()
        void addEnergyPotentials(long *unary_energy_array,
                                 long *binary_energy_array)
        size_t nodeIndex(size_t *node_coords)
        size_t unary_E0_index(size_t node_index)
        size_t unary_E1_index(size_t node_index)
        size_t binary_E00_index(size_t node_index, size_t edge_index)
        size_t binary_E01_index(size_t node_index, size_t edge_index)
        size_t binary_E10_index(size_t node_index, size_t edge_index)
        size_t binary_E11_index(size_t node_index, size_t edge_index)
        void run()
        double runTime()

        vector[int] getCut()

    ctypedef LatticeFlowInterface* solver_ptr

    solver_ptr getLatticeFlowInterface(string solver,
                                       string kernel, size_t* dimensions)

    bint isValidSolver(string name)
    bint isValidKernel(string name)
    size_t kernelDimension(string name)
    vector[string] validKernelNames()
    vector[string] validSolverNames()

cdef extern from "kernels/kernels.hpp" namespace "latticeQBP":
    cdef cppclass Star2d_4:
        pass
    
cdef extern from "tv/tv_solver.hpp" namespace "latticeQBP":

    cdef cppclass FuncArray "latticeQBP::LatticeArray<double, 2>":
        double at(size_t, size_t) nogil except +

    cdef cppclass ImageArray "latticeQBP::LatticeArray<double, 3>":
        double at(size_t, size_t, size_t) nogil except +

    cdef cppclass FuncPathArray "latticeQBP::LatticeArray<double, 3>":
        double at(size_t, size_t, size_t) nogil except +

    ctypedef FuncArray* FuncMapPtr
    ctypedef ImageArray* ImageMapPtr
    ctypedef FuncPathArray* RegPathPtr

    FuncMapPtr _calculate2dTV "latticeQBP::calculate2dTV<latticeQBP::Star2d_24, int64_t>" (
        size_t nx, size_t ny, double *function, double lm) nogil except +

    FuncMapPtr _calculate2dPWReg "latticeQBP::run2dTVProximal<latticeQBP::Star2d_8, int64_t>" (
        size_t nx, size_t ny, double *function, double lm) nogil except +

    ImageMapPtr _calculate2dTVImage "latticeQBP::calculate2dTVImage<latticeQBP::Star2d_24, int64_t>" (
        size_t nx, size_t ny, 
        double *function, double lm,
        double *image, size_t k) nogil except +

    RegPathPtr _calculate2dTV "latticeQBP::calculate2dTV<latticeQBP::Star2d_24, int64_t>" (
        size_t nx, size_t ny, double *function, vector[double]) nogil except +

    
cdef extern from "math.h":
    double exp(double) nogil 

cpdef getKernelDimension(str kernel):
    cdef vector[string] names

    if not isValidKernel(kernel):
        names = validKernelNames()

        raise ValueError(
            "Kernel '%s' not a valid kernel name!  Valid kernels are: %s."
            % (kernel, ', '.join([s for s in names])))

    return kernelDimension(kernel)

cdef class LatticeFlow(object):
    """
    An interface to the lattice flow stuff

    """

    cdef solver_ptr lattice_flow
    cdef size_t n_nodes
    cdef tuple _dimensions
    cdef double numerical_precision_factor
    cdef dict _cache
    cdef ar _edges

    def __init__(self, str kernel, str solver, _dimensions, numerical_precision = 1e-6):

        cdef vector[string] names

        if not isValidSolver(solver):
            names = validSolverNames()

            raise ValueError(
                "'%s' not a valid solver name!  Valid solvers are: %s."
                % (solver, ', '.join(s for s in names)))

        cdef size_t dim = getKernelDimension(kernel)

        self._dimensions = tuple(_dimensions)
        self.n_nodes = np.prod(self._dimensions)

        if len(self._dimensions) != dim:
            raise ValueError("Kernel '%s' is %d-dimensional; a shape of %d dimensions was provided."
                             % (kernel, dim, len(self._dimensions)))

        cdef size_t dimensions[16]
        cdef size_t i, j

        for i in range(dim):
            dimensions[i] = _dimensions[i]

        # Now get the lattice flow set up
        self.lattice_flow = getLatticeFlowInterface(solver, kernel, &dimensions[0])
        self.numerical_precision_factor = 1.0 / numerical_precision

        ############################################################
        # Get the number of edges

        cdef size_t n_edges = self.lattice_flow.nFillingEdges()
        cdef size_t n_dim   = self.lattice_flow.nDimensions()

        cdef ar[int, ndim=2] edges = np.empty( (n_edges, n_dim), dtype='i')

        for i in range(n_edges):
            for j in range(n_dim):
                edges[i,j] = self.lattice_flow.edgeDeltas()[i*n_dim + j]

        self._edges = edges

    @property
    def ndim(self):
        return self.lattice_flow.nDimensions()

    @property
    def shape(self):
        return self._dimensions

    @property
    def size(self):
        return self.n_nodes

    @property
    def nFillingEdges(self):
        return self.lattice_flow.nFillingEdges()

    @property
    def edges(self):
        return self._edges

    def addEnergy(self, ar unary_energy, ar pairwise_energy):

        cdef double nf = self.numerical_precision_factor

        unary_energy = np.ascontiguousarray(np.round(nf * unary_energy), dtype = long)
        pairwise_energy = np.ascontiguousarray(np.round(nf * pairwise_energy), dtype = long)

        cdef long[:, ::1]  unary_a
        cdef long[:, :, ::1] pairwise_a

        ################################################################################
        # Check the unary energy error

        cdef str unary_energy_error = (
            "Unary energy must be either a 2d array with first dimension "
            "matching the number of nodes and second dimension containing E0 and E1, "
            "or (#dimensions + 1) with the first matching the lattice coordinates and "
            "the last containing E0 and E1.")

        if unary_energy.ndim == 2:
            unary_a = unary_energy

            if unary_a.shape[0] != self.n_nodes:
                raise ValueError(unary_energy_error
                                 + (" (# nodes = %d; # unary_array_nodes = %d)"
                                    % (self.n_nodes, unary_a.shape[0])))

            if unary_a.shape[1] != 2:
                raise ValueError(unary_energy_error
                                 + (" (Second dimension has shape %d; 2 required.)"
                                    % (unary_a.shape[1])))


        elif unary_energy.ndim == (len(self._dimensions) + 1):

            for i in xrange(len(self._dimensions)):
                if unary_energy.shape[i] != self._dimensions[i]:
                    raise ValueError(unary_energy_error
                                     + (" (dimension %d = %d; unary_array shape is %d)"
                                        % (i, self._dimensions[i], unary_energy.shape[i])))

            if unary_energy.shape[-1] != 2:
                raise ValueError(unary_energy_error
                                 + (" (Last dimension has shape %d; 2 required.)"
                                    % (unary_energy.shape[-1])))

            unary_a = unary_energy.reshape(-1, 2)

        else:
            raise ValueError(unary_energy_error + " (wrong number of dimensions)")

        ################################################################################
        # Now check the pairwise array

        cdef str pairwise_energy_error = (
            "Pairwise energy must be either: (1) a 3d array with first dimension "
            "matching the number of nodes, second dimension containing the edge direction output, "
            "and third dimension containing E00, E01, E10, and E11 (with E01 being energy of "
            "source node off and dest node on.), or (2) a (#dimensions "
            "+ 2) array with last two dimensions being the same and first dimensions corresponding "
            "to the lattice coordinates. ")

        cdef str bad_edge_error = (
            pairwise_energy_error +
            " (wrong number of edges specified in second-to-last dimension; got %d, needed %d).")

        cdef str bad_last_error = (
            pairwise_energy_error +
            " (wrong number of pairwise energy terms; needed 4, got %d.")


        if pairwise_energy.ndim == 3:
            pairwise_a = pairwise_energy

            if pairwise_a.shape[0] != self.n_nodes:
                raise ValueError(pairwise_energy_error
                                 + (" (# nodes = %d; # nodes in pairwise_energy = %d)"
                                    % (self.n_nodes, unary_a.shape[0])))

            if pairwise_a.shape[1] != self.lattice_flow.nFillingEdges():
                raise ValueError(bad_edge_error
                                 % (pairwise_a.shape[1], self.lattice_flow.nFillingEdges()))

            if pairwise_a.shape[2] != 4:
                raise ValueError(bad_last_error % pairwise_a.shape[2])

        elif pairwise_energy.ndim == (len(self._dimensions) + 2):

            for i in xrange(len(self._dimensions)):
                if pairwise_energy.shape[i] != self._dimensions[i]:
                    raise ValueError(pairwise_energy_error
                                     + (" (dimension %d shape is %d; pairwise_array shape is %d)"
                                        % (i, self._dimensions[i], pairwise_energy.shape[i])))


            if pairwise_energy.shape[-2] != self.lattice_flow.nFillingEdges():
                raise ValueError(bad_edge_error
                                 % (pairwise_energy.shape[-2], self.lattice_flow.nFillingEdges()))

            if pairwise_energy.shape[-1] != 4:
                raise ValueError(bad_last_error % pairwise_energy.shape[-1])

            pairwise_a = pairwise_energy.reshape(-1,self.lattice_flow.nFillingEdges(), 4)

        else:
            raise ValueError(pairwise_energy_error + " (Array has wrong number of dimensions.)")

        ################################################################################
        # And now we can just call it

        self.lattice_flow.addEnergyPotentials(&unary_a[0,0], &pairwise_a[0,0,0])

    def run(self):
        """
        Runs the specified solver.
        """

        self.lattice_flow.run()


    def runTime(self):
        """
        Returns the runing time of the solver (the last call to run) in seconds.
        """

        return self.lattice_flow.runTime();


    def getCut(self):
        """
        Gets the value of the cut as an integer array with dimensions
        matching the lattice.
        """

        cdef vector[int] cut = self.lattice_flow.getCut()

        assert cut.size() == self.n_nodes

        cdef ar[int, mode='c'] cut_a = np.empty(cut.size(), dtype = 'i')

        cdef size_t i

        for i in range(cut.size()):
            cut_a[i] = cut[i]

        return cut_a.reshape(*self._dimensions)


def _runBenchmark(str image_file, str kernel, str solver, bint check_result = False):

    from matplotlib.pylab import imread

    Xo = imread(image_file)

    if not Xo.size:
        raise IOError("Error loading image %s." % image_file)

    if len(Xo.shape) != 3:
        raise IOError("Error loading image %s; need RGB images." % image_file)

    cdef double[:,:,:] X = np.asarray(Xo[:,:,:3], dtype='d')

    cdef size_t n_nodes

    cdef size_t dimension = getKernelDimension(kernel)

    if dimension == 2:
        n_nodes = X.shape[0] * X.shape[1]
        dimensions = (X.shape[0], X.shape[1])
    elif dimension == 3:
        n_nodes = X.shape[0] * X.shape[1] * X.shape[1]
        dimensions = (X.shape[0], X.shape[1], X.shape[1])
    else:
        raise ValueError("Benchmark for dimension %d not implemented (kernel=%s, solver=%s)"
                         % (dimension, kernel, solver))

    cdef LatticeFlow lf = LatticeFlow(kernel, solver, dimensions)

    # Build the lists
    cdef int[:,:] delta_list = lf.edges

    cdef double[:] scale_factors = np.empty(delta_list.shape[0])

    cdef size_t i, j, xi, yi, zi, ei, xj, yj, zj, node_idx

    for i in range(scale_factors.shape[0]):
        scale_factors[i] = 0

        for j in range(0, delta_list.shape[1]):
            scale_factors[i] += delta_list[i, j]*delta_list[i, j]

        scale_factors[i] = 1.0 / scale_factors[i];

    cdef ar[double, ndim=2,mode='c'] unary_energies = np.zeros( (n_nodes, 2) )
    cdef ar[double, ndim=3,mode='c'] binary_energies = np.zeros( (n_nodes, delta_list.shape[0], 4) )
    cdef double d, x1, x2, y1, y2, z1, z2

    cdef bint mode_3d = (dimension == 3)

    print "Calculating potentials."

    for xi in range(X.shape[0]):

        for yi in range(X.shape[1]):

            for zi in range(X.shape[1] if mode_3d else 1):

                node_idx = (xi*X.shape[1]*X.shape[1] + yi*X.shape[1] + zi) if mode_3d else (xi*X.shape[1] + yi)
                # unary_energies[node_idx, 0] = 64 - ((32 * <double>(xi +1 + yi)* (xi + 3) * (4*yi + 1) + 50 * zi + 10*xi) % 93)
                # continue

                if not mode_3d:
                    x1, y1, z1 = X[xi,yi,0], X[xi,yi,1], X[xi,yi,2]
                else:
                    x1 = (X[xi, yi, 0] + X[xi, zi, 0]) / 2
                    y1 = (X[xi, yi, 1] + X[xi, zi, 1]) / 2
                    z1 = (X[xi, yi, 2] + X[xi, zi, 2]) / 2

                if x1 == y1 == z1 == 0:
                    unary_energies[node_idx, 0] = 50000
                elif x1 == y1 == z1 == 1:
                    unary_energies[node_idx, 1] = 50000

                for ei in range(delta_list.shape[0]):

                    xj = xi + delta_list[ei, 0]
                    yj = yi + delta_list[ei, 1]
                    zj = zi + (delta_list[ei, 2] if mode_3d else 0)

                    if xj >= X.shape[0] or yj >= X.shape[1] or (mode_3d and zj >= X.shape[1]):
                        continue

                    if not mode_3d:
                        x2, y2, z2 = X[xj,yj,0], X[xj,yj,1], X[xj,yj,2]
                    else:
                        x2 = (X[xj, yj, 0] + X[xj, zj, 0]) / 2
                        y2 = (X[xj, yj, 1] + X[xj, zj, 1]) / 2
                        z2 = (X[xj, yj, 2] + X[xj, zj, 2]) / 2

                    d = exp(-(((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)))*scale_factors[ei]

                    # Set the breaking penalties
                    binary_energies[node_idx, ei, 1] = d
                    binary_energies[node_idx, ei, 2] = d

    min_d = binary_energies[:,:,1].min()
    binary_energies[:,:,1:2] -= min_d

    lf.addEnergy(unary_energies, binary_energies)

    lf.run()

    cut = lf.getCut().ravel()
    runtime = lf.runTime()

    print ("Final result: %d / %d turned on, solution hash = %s, ran in %f seconds."
           % (cut.sum(), cut.shape[0], md5(dumps(cut)).hexdigest()[:8], runtime))

    if check_result:
        print "Checking solution:"

        lf = LatticeFlow(kernel, "BKGCEnergyMinimizer", dimensions)
        lf.addEnergy(unary_energies, binary_energies)

        lf.run()

        print "BKGC Solver ran in %f seconds." % lf.runTime()

        cut2 = lf.getCut().ravel()

        diff = (cut != cut2).sum()

        if diff.sum() != 0:
            print ("Solution different: %d falsely 0, %d falsely 1"
                   % (((cut == 0) & (cut2 == 1)).sum(), ((cut == 1) & (cut2 == 0)).sum()))
        else:
            print "Solutions exactly match."

    return lf.runTime()

def calculate2dTV(ar Xo, double flt):

    cdef ar[double, ndim=2, mode='c'] X = np.ascontiguousarray(Xo, dtype='d')

    cdef size_t nx = X.shape[0]
    cdef size_t ny = X.shape[1]

    cdef FuncMapPtr Rv = _calculate2dTV(nx, ny, &X[0,0], flt)

    cdef ar[double, ndim=2, mode='c'] R = np.empty( (X.shape[0], X.shape[1]) )

    cdef size_t i, j

    for xi in range(nx):
        for yi in range(ny):
            R[xi,yi] = Rv.at(xi, yi)
        
    return R

def calculate2dTVImage(ar Xo, double flt):

    cdef ar[double, ndim=3, mode='c'] X = np.ascontiguousarray(Xo, dtype='d')

    cdef ar[double, ndim=2, mode='c'] X_r = (X.mean(axis=2) / X.max())

    cdef size_t nx = X.shape[0]
    cdef size_t ny = X.shape[1]
    cdef size_t nz = X.shape[2]

    cdef ImageMapPtr Rv = _calculate2dTVImage(nx, ny, &X_r[0,0], flt, &X[0,0,0], nz)

    cdef ar[double, ndim=3, mode='c'] R = np.empty( (X.shape[0], X.shape[1], X.shape[2]) )

    cdef size_t i, j, k

    for xi in range(nx):
        for yi in range(ny):
            for zi in range(nz):
                R[xi,yi,zi] = Rv.at(xi, yi,zi)
        
    return R

    
def calculate2dTVPath(ar Xo, ar[double] lma):

    cdef size_t i, j

    cdef ar[double, ndim=2, mode='c'] X = np.ascontiguousarray(Xo, dtype='d')
    
    cdef size_t nx = X.shape[0]
    cdef size_t ny = X.shape[1]

    cdef size_t n_lm = lma.shape[0]

    cdef vector[double] lm_values = vector[double](n_lm)

    for i in range(n_lm):
        lm_values[i] = lma[i]

    cdef RegPathPtr Rv = _calculate2dTV(nx, ny, &X[0,0], lm_values)

    cdef ar[double, ndim=3, mode='c'] R = np.empty( (n_lm, X.shape[0], X.shape[1]) )

    cdef size_t li, xi, yi

    for li in range(n_lm):
        for xi in range(nx):
            for yi in range(ny):
                R[li,xi,yi] = Rv.at(li, xi, yi)
        
    return R


@cython.boundscheck(False)
cdef double calculate_g(ar[double, ndim=2, mode="c"] X, double reg_p):
    cdef int i, j, ii, jj

    cdef double r = 0, gc
    
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):

            # ii, jj, gc = i + -3, j + -1, 0.0366545619033433
            # if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
            #     r += gc*fabs(X[i,j] - X[ii,jj])

            # ii, jj, gc = i + -3, j + 1, 0.0366545619033433
            # if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
            #     r += gc*fabs(X[i,j] - X[ii,jj])

            # ii, jj, gc = i + -2, j + -1, 0.0518373785665523
            # if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
            #     r += gc*fabs(X[i,j] - X[ii,jj])

            # ii, jj, gc = i + -2, j + 1, 0.0518373785665523
            # if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
            #     r += gc*fabs(X[i,j] - X[ii,jj])

            # ii, jj, gc = i + -1, j + -3, 0.0366545619033433
            # if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
            #     r += gc*fabs(X[i,j] - X[ii,jj])

            # ii, jj, gc = i + -1, j + -2, 0.0518373785665523
            # if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
            #     r += gc*fabs(X[i,j] - X[ii,jj])

            # ii, jj, gc = i + -1, j + -1, 0.1137559994321984
            # if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
            #     r += gc*fabs(X[i,j] - X[ii,jj])

            # ii, jj, gc = i + -1, j + 0, 0.1608752771983211
            # if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
            #     r += gc*fabs(X[i,j] - X[ii,jj])

            # ii, jj, gc = i + -1, j + 1, 0.1137559994321984
            # if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
            #     r += gc*fabs(X[i,j] - X[ii,jj])

            # ii, jj, gc = i + -1, j + 2, 0.0518373785665523
            # if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
            #     r += gc*fabs(X[i,j] - X[ii,jj])

            # ii, jj, gc = i + -1, j + 3, 0.0366545619033433
            # if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
            #     r += gc*fabs(X[i,j] - X[ii,jj])

            # ii, jj, gc = i + 0, j + -1, 0.1608752771983211
            # if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
            #     r += gc*fabs(X[i,j] - X[ii,jj])

            # ii, jj, gc = i + 0, j + 1, 0.1608752771983211
            # if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
            #     r += gc*fabs(X[i,j] - X[ii,jj])

            # ii, jj, gc = i + 1, j + -3, 0.0366545619033433
            # if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
            #     r += gc*fabs(X[i,j] - X[ii,jj])

            # ii, jj, gc = i + 1, j + -2, 0.0518373785665523
            # if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
            #     r += gc*fabs(X[i,j] - X[ii,jj])

            # ii, jj, gc = i + 1, j + -1, 0.1137559994321984
            # if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
            #     r += gc*fabs(X[i,j] - X[ii,jj])

            # ii, jj, gc = i + 1, j + 0, 0.1608752771983211
            # if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
            #     r += gc*fabs(X[i,j] - X[ii,jj])

            # ii, jj, gc = i + 1, j + 1, 0.1137559994321984
            # if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
            #     r += gc*fabs(X[i,j] - X[ii,jj])

            # ii, jj, gc = i + 1, j + 2, 0.0518373785665523
            # if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
            #     r += gc*fabs(X[i,j] - X[ii,jj])

            # ii, jj, gc = i + 1, j + 3, 0.0366545619033433
            # if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
            #     r += gc*fabs(X[i,j] - X[ii,jj])

            # ii, jj, gc = i + 2, j + -1, 0.0518373785665523
            # if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
            #     r += gc*fabs(X[i,j] - X[ii,jj])

            # ii, jj, gc = i + 2, j + 1, 0.0518373785665523
            # if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
            #     r += gc*fabs(X[i,j] - X[ii,jj])

            # ii, jj, gc = i + 3, j + -1, 0.0366545619033433
            # if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
            #     r += gc*fabs(X[i,j] - X[ii,jj])

            # ii, jj, gc = i + 3, j + 1, 0.0366545619033433
            # if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
            #     r += gc*fabs(X[i,j] - X[ii,jj])

            # ii, jj, gc = i + 3, j + -1, 0.0366545619033433
            # if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
            #     r += gc*fabs(X[i,j] - X[ii,jj])

            # ii, jj, gc = i + 3, j + 1, 0.0366545619033433
            # if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
            #     r += gc*fabs(X[i,j] - X[ii,jj])

            # ################################################################################

            ii, jj, gc = i + 0, j - 1, 0.3926990816987241,
            if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
                r += gc*fabs(X[i,j] - X[ii,jj])

            ii, jj, gc = i + 1, j - 1, 0.2776801836348979
            if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
                r += gc*fabs(X[i,j] - X[ii,jj])

            ii, jj, gc = i + 1, j + 0, 0.3926990816987241
            if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
                r += gc*fabs(X[i,j] - X[ii,jj])

            ii, jj, gc = i + 1, j + 1, 0.2776801836348979
            if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
                r += gc*fabs(X[i,j] - X[ii,jj])

            ################################################################################


            # ii, jj, gc = i + 0, j + 1, 0.7853981633974483
            # if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
            #     r += gc*fabs(X[i,j] - X[ii,jj])

            # ii, jj, gc = i + 1, j + 0, 0.7853981633974483
            # if 0 <= ii < X.shape[0] and 0 <= jj < X.shape[1]:
            #     r += gc*fabs(X[i,j] - X[ii,jj])

    return  2*r * reg_p


cdef QL(ar[double, ndim=2, mode="c"] X,
               ar[double, ndim=2, mode="c"] Y,
               double value_Y,
               ar[double, ndim=2, mode="c"] grad_Y,
               double L,
               double value_g):

    XYdiff = X - Y
    
    return (value_Y 
            + ((XYdiff)*grad_Y).sum()
            + 0.5*L * (la.norm(XYdiff)**2)
            + value_g)

def check_shit(Yo):

    cdef ar[double, ndim=2, mode="c"] Y = Yo

    print "Y = "
    print Y

    cdef size_t nx = Y.shape[0]
    cdef size_t ny = Y.shape[1]

    cdef FuncMapPtr Rv = _calculate2dPWReg(nx, ny, &Y[0,0], 0)

    cdef size_t xi, yi
    cdef size_t i, j

    for xi in range(nx):
        for yi in range(ny):
            assert fabs(Rv.at(xi, yi) - Y[xi, yi]) < 1e-5

    cdef double lm = 0.001

    Rv = _calculate2dPWReg(nx, ny, &Y[0,0], lm)

    cdef ar[double, ndim=2, mode='c'] R = np.empty( (Y.shape[0], Y.shape[1]) )
    
    for xi in range(nx):
        for yi in range(ny):
            R[xi,yi] = Rv.at(xi, yi)

    # print "R = "
    # print R

    # print " "

    cdef double v1, v2, v3, cv, gv1, gv2, gv3, fv1, fv2, fv3

    cdef double eps = .1
            
    # Make sure it's found the optimal there
    fv1 = ((Y - R)**2).sum()
    gv1 = calculate_g(R, lm)
    v1 = gv1 + fv1

    assert gv1 != 0

    cdef ar[double, ndim=2] R2 = R.copy()

    is_okay = True

    Lset = set(list(np.linspace(0,16,50)) + list(np.linspace(-16,0,50)))

    for xi in range(nx):
        for yi in range(ny):

            R2[xi,yi] += eps

            gv2 = calculate_g(R2, lm)

            assert gv2 != 0
            fv2 = ((Y - R2)**2).sum()

            v2 = fv2 + gv2

            R2[xi, yi] -= 2*eps

            gv3 = calculate_g(R2, lm)
            fv3 = ((Y - R2)**2).sum()

            v3 = fv3 + gv3

            if v2 < v1 or v3 < v1:
                L = [c for c in np.linspace(0,16,500)
                     if ((fv1 + (2.0**c) * gv1) <= (fv2 + (2.0**c) * gv2) 
                         and (fv1 + (2.0**c) * gv1) <= (fv3 + (2.0**c) * gv3))]

                L += [c for c in np.linspace(-16,0,500)
                      if ((fv1 + (2.0**c) * gv1) <= (fv2 + (2.0**c) * gv2) 
                         and (fv1 + (2.0**c) * gv1) <= (fv3 + (2.0**c) * gv3))]

                Lset = Lset & set(L)
                
                print "gv2=%f, gv1=%f, gv3=%f" % (gv2,gv1,gv3)
                print "fv2=%f, fv1=%f, fv3=%f" % (fv2,fv1,fv3)
                print "v2=%f,  v1=%f,  v3=%f" % (v2,v1,v3)
                print " " 

                is_okay = False

            R2[xi, yi] += eps


    if not is_okay:
        print "Lset = ", "\n".join([("%f" % v) for v in sorted(Lset)])
        assert is_okay
    



cdef pL(ar[double, ndim=2, mode="c"] Y,
        ar[double, ndim=2, mode="c"] grad_Y,
        double L, double reg_p, f):

    cdef ar[double, ndim=2, mode="c"] Yt = Y - (1.0 / L) * grad_Y

    # check_shit(Yt)


    cdef size_t nx = Yt.shape[0]
    cdef size_t ny = Yt.shape[1]

    # print "Yt = ",
    # print Yt

    cdef FuncMapPtr Rv = _calculate2dPWReg(nx, ny, &Yt[0,0], 2 * reg_p / L)

    cdef ar[double, ndim=2, mode='c'] R = np.empty( (Y.shape[0], Y.shape[1]) )

    cdef size_t xi, yi
    cdef size_t i, j

    for xi in range(nx):
        for yi in range(ny):
            R[xi,yi] = Rv.at(xi, yi)

    # print "R = ",
    # print R
    
    return R
    
    


@cython.boundscheck(False)
cpdef _run2DFista(f, ar[double, ndim=2, mode="c"] initial_X,
                  double reg_p,
                  size_t n_iterations, size_t sample_interval, double eta = 1.3,
                  double diff_eps = 1e-6):

    assert eta > 1

    cdef ar[double, ndim=2, mode="c"] X = initial_X.copy()

    cdef double L = 1

    cdef size_t iteration

    cdef double t = 1, t_new = 1

    cdef double value_f_Xtr, value_g_Xtr, value_F_Xtr, value_QL
    
    cdef list solution_track = []

    #print "X.shape = ", X.shape[0], X.shape[1]

    Y = X.copy()

    for iteration in range(n_iterations):

        # print "Y.shape = ", Y.shape[0], Y.shape[1]
        
        value_Y, grad_Y = f(Y)

        # value_Y *= reg_p
        # grad_Y *= reg_p

        assert grad_Y.shape == Y.shape

        if sample_interval != 0 and (iteration % sample_interval) == 0:
            value_F_Y = value_Y + calculate_g(Y, reg_p)
            solution_track.append( (value_F_Y, Y.copy()) )

            print "Iteration", iteration, "; value = ", value_F_Y

            if len(solution_track) >= 5 and np.std([v for v, s in solution_track[-5:]]) <= diff_eps:
                break
            
        # print "; ".join(["%1.5f" % v for v in (Y.ravel())])

        # print "Y = "
        # print Y

        while True:
            
            Xtr = pL(Y, grad_Y, L, reg_p, f)

            # print "Xtr = " #"; ".join(["%1.5f" % v for v in (Y.ravel())])
            # print Xtr

            value_f_Xtr, grad_f_Xtr = f(Xtr)
            
            value_g_Xtr = calculate_g(Xtr, reg_p)
            value_F_Xtr = value_f_Xtr + value_g_Xtr
            
            value_QL = QL(Xtr, Y, value_Y, grad_Y, L, value_g_Xtr)

            # print " >>>> value_QL, value_F_Xtr = ", value_QL, value_F_Xtr

            if value_F_Xtr > value_QL:
                L *= eta
                print ">>> L = ", L 
            else:
                break

        t_new = 0.5*(1 + fsqrt(1 + 4*t*t))

        Y = Xtr + ((t - 1) / t_new ) * (Xtr - Y)

        t = t_new

    return Y, solution_track


def regularizedRegression(A, ar y, double reg_parameter,
                          size_t n_iterations = 1000, size_t sample_interval = 10,
                          double eta = 1.2, double diff_eps = 0):
    
    """
    Runs that model.  Blah!

    """

    assert A.ndim == 3

    Ars = A.reshape(-1, A.shape[-1]).T

    yA = dot(y, Ars).reshape(1,-1)
    
    def f(X):

        assert X.shape == A.shape[:-1]

        Xrs = X.reshape(-1,1)

        assert yA.size == Xrs.shape[0]
        assert Ars.shape[1] == Xrs.shape[0], \
               ("A.shape[1] = %d != %d = Xrs.shape[0]"
                % (A.shape[1], Xrs.shape[0]))

        xtAtA = dot(dot(Xrs.T, Ars.T), Ars)

        # print "xtAtA.shape", xtAtA.shape
        # print "Xrs.shape", Xrs.shape
        # print "yA.shape", yA.shape

        value = dot(xtAtA, Xrs) - 2 * dot(yA, Xrs)

        assert value.size == 1

        grad = 2 * (xtAtA - yA)

        # print "Grad = "
        # print grad.reshape(X.shape)

        return value[0,0], grad.reshape(X.shape)

    starting_X = (dot(la.pinv(Ars), y.reshape(-1,1))).reshape(A.shape[:2])

    #dot(la.pinv(Ars), y).reshape(A.shape[:-1])

    print "Running Fista"

    return _run2DFista(f, starting_X, reg_parameter,
                       n_iterations, sample_interval, eta, diff_eps)
