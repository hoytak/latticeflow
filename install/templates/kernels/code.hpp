struct %(kernel_name)s : public KernelBase<%(dimension)d, %(n_edges)d, %(is_geocut_applicable)d> {
  %(kernel_name)s() 
    : KernelBase<%(dimension)d, %(n_edges)d, %(is_geocut_applicable)d>
    (%(edge_list)s, 
     %(geocut_edge_weights)s
    )
  {}
};
