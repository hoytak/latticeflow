#ifndef _SPLIT_HEAP_H_
#define _SPLIT_HEAP_H_

template <typename T, typename ComparisonFunctor>
class __ReverseComparisonFunctor {
public:
  bool operator()(const T& t1, const T& t2) const {
    return ComparisonFunctor(t2, t1);
  }
};

template <typename T, typename ComparisonFunctor>
class SplitHeap {
public:
  
  inline const T& u_bottom() const {
    return u_heap.front();
  }

  inline const T& l_top() const {
    return l_heap.front();
  }
  
  inline void moveFromLtoU() {
    pop_heap(l_heap.begin(), l_heap.end(), __ReverseComparisonFunctor());
    u_heap.push_back(l_heap.back());
    l_heap.pop_back();
    push_heap(u_heap.begin(), u_heap.end(), ComparisonFunctor);
  }

  inline void moveFromUtoL() {
    pop_heap(u_heap.begin(), u_heap.end(), ComparisonFunctor);
    l_heap.push_back(u_heap.back());
    u_heap.pop_back();
    push_heap(l_heap.begin(), l_heap.end(), __ReverseComparisonFunctor());
  }
  
  inline void addToL(const T& t) {
    l_heap.push_back(t);
    push_heap(l_heap.begin(), l_heap.end(), __ReverseComparisonFunctor());
  }

  inline void addToU(const T& t) {
    u_heap.push_back(t);
    push_heap(u_heap.begin(), u_heap.end(), ComparisonFunctor);
  }

  inline void popFromL() {
    pop_heap(l_heap.begin(), l_heap.end(), __ReverseComparisonFunctor());
    l_heap.pop_back();
  }

  inline void popFromU() {
    pop_heap(u_heap.begin(), u_heap.end(), ComparisonFunctor);
    u_heap.pop_back();
  }

private:
  vector<T> l_heap, h_heap;
};

#endif /* _CRITICAL_POINT_HEAP_H_ */
