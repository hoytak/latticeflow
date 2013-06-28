struct %(kernel_name)s : public KernelBase<%(dimension)d, %(n_edges)d> {
%(kernel_name)s() : KernelBase<%(dimension)d, %(n_edges)d>(
  %(edge_list)s)
  {}
};
