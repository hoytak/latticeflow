// Only meant to be included from within energy.hpp

#ifndef _ENERGY_PRIVATE_MEMBER_
#error "File nodes.hpp cannot be included directly."
#else

template <typename Kernel, int only_forward_edges>
class NetworkFlowNode {
  Node()
    : state(DEBUG_MODE ? -1 : 0)
    , height(0)
    , level_index(0)
    , reduction(0)
    , alpha(0)
  {}

  int state;

  // The height for the push-relabel thing
  unsigned int height;

  unsigned int level_index;
  ftype reduction;
  
  constexpr int nEdges() {
    return only_forward_edges ? kernel::kernelPositiveSize() : kernel::kernelSize(); 
  }

  Array<dtype, nEdges()> alpha;

  //////////////////////////////////////////////////
  // Helper methods that tie in other things...
  bool on() const { 
    assert(state == 1 || state == -1);
    return (state == 1);
  }

  // Note that excess, in this case, is defined in relation to the
  // push-relabel algorithm; so excess-flow corresponds to deficit in
  // the reduction.  

  template <int partition=0> inline ftype excess() const {
    assert_equal(partition, on());
    return ((partition == 1) ? 1 : -1) * (reduction);
  }

  template <int partition=0> inline ftype sinkFlow() const {
    return -excess<partition>();
  }

  template <int partition=0> inline ftype pullCapacity(int dim) const {
    assert_leq(alpha[dim], -edges[dim]);
    assert_geq(alpha[dim],  edges[dim]);

    // Actually the amount that can be pushed to the other node, so in
    // the 1-partition it's 

    return -edges[dim] + ( (partition == 1) ? 1 : -1) * alpha[dim];
  }

  template <int partition> inline ftype pushCapacity(int dim) const {
    return -edges[dim] - ((partition == 1) ? 1 : -1) * alpha[dim];
  }

  template <int partition> inline bool pushSaturated(int dim) const {
    return (pushCapacity<partition>(dim) == 0);
  }

  template <int partition> inline bool pullSaturated(int dim) const {
    return (pullCapacity<partition>(dim) == 0);
  }

  inline bool isConnected(int dim) const {
    return edges[dim] != 0;
  }

  static dtype pinNodeConstant() {
    return (1 << (8*sizeof(dtype) - 3));
  }

  inline void pin() {
    assert(!isPinned());
    gain -= pinNodeConstant();
    assert(isPinned());
  }
  
  inline void unpin() {
    assert(isPinned());
    gain += pinNodeConstant();
    assert(!isPinned());
  }

  inline bool isPinned() const {
    return (gain <= -(pinNodeConstant() / 2));
  }
};

template <int partition> 
inline ftype pushCapacity(const node_cptr& src, const node_cptr& dest, const DirDim& dd) const {
  
  assert(dest == lattice.neighbor(src, dd));

  ftype v;

  if(dd.direction() == 1) {
    v = src->template pushCapacity<partition>(dd.dimension());
  } else {
    v = dest->template pullCapacity<partition>(dd.dimension());
  }
  
  assert_geq(v, 0);
  
  // cout << "Push capacity from node " << (src - lattice.begin()) 
  //      << " to node " << (dest - lattice.begin()) << " is " 
  //      << v << ". " << endl;

  return v;
}

template <int partition>
inline bool isSaturated(const node_cptr& src, const node_cptr& dest, const DirDim& dir_dim) const {
  return (pushCapacity<partition>(src, dest, dir_dim) == 0);
}

template <int partition> 
inline ftype pushCapacity(const node_cptr& src, const node_cptr& dest, int dim, int dir) const {
  
  assert(dest == lattice.neighbor(src, dim, dir));

  ftype v;

  if(dir == 1) {
    v = src->template pushCapacity<partition>(dim);
  } else {
    v = dest->template pullCapacity<partition>(dim);
  }
  
  assert_geq(v, 0);
  
  return v;
}

inline bool isConnected(const node_ptr& src, const node_ptr& dest, const DirDim& dd) {
  int direction = dd.direction();
  int dim = dd.dimension();
  
  assert(dest == lattice.neighbor(src, dim, direction));

  return ((direction > 0) ? src : dest)->isConnected(dim);
}

inline bool isConnected(const node_ptr& src, const node_ptr& dest, int index) {
  int direction = lattice.indexDirection(index);
  int dim = lattice.indexDimension(index);
  
  assert(dest == lattice.neighbor(src, dim, direction));

  return ((direction > 0) ? src : dest)->isConnected(dim);
}

template <int partition> 
inline void pushExcess(const node_ptr& src, const node_ptr& dest, const DirDim& dd, ftype amount) {

  int direction = dd.direction();
  int dim = dd.dimension();

#ifndef NDEBUG

  assert(direction == -1 || direction == 1);
  assert(dest == lattice.neighbor(src, dim, direction));
  assert_equal(src->on(), partition);
  assert_equal(dest->on(), partition);
  assert_geq(dim, 0);
  assert_leq(amount, pushCapacity<partition>(src, dest, dd));
  assert_leq(-amount, pushCapacity<partition>(dest, src, dd.reversed()));
  assert(direction == 1 || direction == -1);
  
  _debugVerifyNodeReduction(src);
  _debugVerifyNodeReduction(dest);

  ftype src_pr_excess = src->template excess<partition>();
  ftype dest_pr_excess = dest->template excess<partition>();

#endif

  src->reduction  -= ((partition == 1) ? 1 : -1) * amount;
  dest->reduction += ((partition == 1) ? 1 : -1) * amount;

  ((direction > 0) ? src : dest)
    ->alpha[dim] += direction * ((partition == 1) ? 1 : -1) * amount;
  
#ifndef NDEBUG
  _debugVerifyNodeReduction(src);
  _debugVerifyNodeReduction(dest);

  ftype src_af_excess = src->template excess<partition>();
  ftype dest_af_excess = dest->template excess<partition>();
  
  assert_equal(src_pr_excess - amount, src_af_excess);
  assert_equal(dest_pr_excess + amount, dest_af_excess);
#endif
  
}

////////////////////////////////////////////////////////////////////////////////
// inter-partition work

template <int start_state>
inline void flipNode(node_ptr node) {
      
  // Go through and remove all the reductions from neighbors that
  // are in the same state as we are; then impose the new
  // reductions on them.

  // Also need to update the base reductions of these other nodes
  // with the influence difference of this new state

#ifndef NDEBUG
  assert(node >= lattice.begin() && node < lattice.end());
  assert(!node->isPinned());

  _debugVerifyNodeConsistency(node);
  _debugVerifyNodeReduction(node);

  for(size_t i = 0; i < n_dimensions; ++i) {
    _debugVerifyNodeConsistency(lattice.neighbor(node, i, -1));
    _debugVerifyNodeConsistency(lattice.neighbor(node, i, 1));
    _debugVerifyNodeReduction(lattice.neighbor(node, i, -1));
    _debugVerifyNodeReduction(lattice.neighbor(node, i, 1));
  }
#endif
      
  // cout << "Flipping Node WR " << (node - lattice.begin()) << endl;

  const int s1 = (start_state <= 1) 
    ? ((start_state == 0) ? -1 : 1)
    : node->state;

  for(size_t i = 0; i < n_dimensions; ++i) {
    
    const dtype node_edge = node->edges[i];
    const dtype node_alpha = node->alpha[i];

    {
      const int dir = 1;
      const node_ptr nn = lattice.neighbor(node, i, dir);
      const int s2 = nn->state;

      const dtype r = dir*s1*s2*node_alpha + s1*node_edge;

      node->reduction += r;
      nn->reduction   -= r;
      nn->gain        -= s1 * s2 * node_edge;
    }

    {
      const int dir = -1;
      const node_ptr nn = lattice.neighbor(node, i, dir);
      const int s2 = nn->state;
      const dtype nn_edge = nn->edges[i];
      const dtype nn_alpha = nn->alpha[i];

      const dtype r    = dir*s1*s2*nn_alpha + s1*nn_edge;

      node->reduction += r;
      nn->reduction   -= r;
      nn->gain        -= s1 * s2 * nn_edge;
    }
  }

  node->gain *= -1;
  current_value += node->gain;

  node->state = -s1;

#ifndef NDEBUG
  assert(!node->isPinned());

  _debugVerifyNodeConsistency(node);
  _debugVerifyNodeReduction(node);

  for(size_t i = 0; i < n_dimensions; ++i) {
    _debugVerifyNodeConsistency(lattice.neighbor(node, i, -1));
    _debugVerifyNodeConsistency(lattice.neighbor(node, i, 1));
    _debugVerifyNodeReduction(lattice.neighbor(node, i, -1));
    _debugVerifyNodeReduction(lattice.neighbor(node, i, 1));
  }
#endif
}


#endif
