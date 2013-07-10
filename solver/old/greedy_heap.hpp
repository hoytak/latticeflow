#ifndef _ENERGY_PRIVATE_MEMBER_
#error "File greedy_heap.hpp cannot be included directly."
#else

class GreedyHeap {
public:

  inline bool top(node_ptr& n, bool must_have_something = false) {
    while(true) {

      if(DEBUG_MODE) {
	if(must_have_something)
	  assert(!queue.empty());
      }

      if(!must_have_something && unlikely(queue.empty()))
	return true;

      const GHNode& ghn = queue_top();
      n = ghn.node;
      dtype gain = ghn.gain;
      
      if(likely(gain == n->gain))
	return false;

      queue_pop();
    }
  }

  inline void pop() {
    queue_pop();
  }

  inline void push(node_ptr n) {
    assert(!n->isPinned());
    queue_push(n);
  }

  inline void clear() {
    queue.clear();
  }

  inline size_t size() const {
    return queue.size();
  }

  inline bool empty() const {
    return queue.empty();
  }
  
private:
  
  struct GHNode {
    GHNode(node_ptr _n) 
      : gain(_n->gain), node(_n)
    {}

    bool operator<(const GHNode& ghn) const {
      return gain < ghn.gain;
    }

    dtype gain;
    node_ptr node;
  };
  
  const GHNode& queue_top() {
    return queue.front();
  }

  void queue_pop() {
    pop_heap(queue.begin(), queue.end());
    queue.pop_back();
  }

  void queue_push(node_ptr n) {
    queue.push_back(GHNode(n));
    push_heap(queue.begin(), queue.end());
  }

  vector<GHNode> queue;
};

#endif
