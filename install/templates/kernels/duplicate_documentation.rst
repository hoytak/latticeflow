
.. _%(kernel_name)s:

%(kernel_header_name)s
----------------------------------------------------------------------

%(description)s


Note: This kernel is identical to :ref:`%(reference_kernel_name)s`.

Properties:

:C++ Name: %(kernel_name)s

:Dimension: %(dimension)d

:Radius: %(radius)s

:Edges: %(n_edges)s


%(edge_graph_header)s:

.. graphviz::

   digraph graph_%(kernel_name)s {
     graph [center, layout=neato, maxiter=0, bgcolor="#C0C0C0", sep=0.2, splines=true]

     edge [dir=none]

     node [width=0.15 height=0.15 label=""] {
%(graphviz_nodes)s
     }

%(graphviz_grid_edges)s

%(graphviz_edges)s
   }
