static solver_ptr new__%(solver_name)s__%(reference_kernel_name)s(size_t *indices) {
  return solver_ptr(
    new LatticeFlowDirectInterface<%(solver_name)s<%(dimension)d, %(reference_kernel_name)s> >
      (Array<size_t, %(dimension)d>(indices))
    );
}
