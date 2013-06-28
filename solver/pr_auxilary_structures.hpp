#ifndef _PR_AUXILARY_STRUCTURES_H_
#define _PR_AUXILARY_STRUCTURES_H_

#include <cstdint>

template <typename level_index_type, int kernel_size, int use_compact_form> 
class PossibleNeighbor {};

template <typename level_index_type, int kernel_size> 
class PossibleNeighbor<level_index_type, kernel_size, 0> {
public:
  PossibleNeighbor() : _height(0), _index(0) {}

  inline bool operator<(const PossibleNeighbor& nq) const {
    return _height < nq._height;
  }
    
  inline void set(level_index_type h, int idx) {
    _height = h;
    _index = idx;
  }

  inline level_index_type height() const {
    return _height;
  }

  inline level_index_type index() const {
    return _index;
  }

  inline void incrementHeight() {
    ++_height;
  }

  inline bool nonzeroHeightEqualsOne() const {
    return _height == 1;
  }

  inline bool heightIsGreaterThan(level_index_type h) const {
    return _height > h;
  }

private:
  level_index_type _height;
  unsigned int _index;
};

template <typename level_index_type, int kernel_size> 
class PossibleNeighbor<level_index_type, kernel_size, 1> {
public:
  static inline int shiftAmount() {
    unsigned int i = 0, k = kernel_size;
    while(k != 0) {
      k >>= 1;
      ++i;
    }
    return i;
  }

  static inline uint64_t _dataValue(level_index_type h, int idx) {
    return (uint64_t(h) << shiftAmount()) + idx;
  }

  PossibleNeighbor() : data(0) {}

  PossibleNeighbor(level_index_type h, int idx) 
    : data(_dataValue(h, index)) 
  {
    assert_equal(h, height());
    assert_equal(idx, index());
  }

  inline PossibleNeighbor(const PossibleNeighbor& n)
    : data(n.data)
  {}

  inline const PossibleNeighbor& operator=(const PossibleNeighbor& n) {
    data = n.data;
    return *this;
  }

  inline void set(level_index_type h, int idx) {
    data = _dataValue(h, idx);

    assert_equal(h, height());
    assert_equal(idx, index());
  }

  inline bool operator<(const PossibleNeighbor& nq) const {
    return data < nq.data;
  }

  inline void incrementHeight() {
    data += (1 << shiftAmount());
  }
    
  inline bool nonzeroHeightEqualsOne() const {
    assert_gt(height(), 0);
    return !(data >> (1 + shiftAmount()));
  }

  inline bool heightIsGreaterThan(level_index_type h) const {
    // Easier for compiler to move computations out of inner loops
    return (data >= ( (h + 1) << shiftAmount())); 
  }

  inline level_index_type height() const { 
    return level_index_type(data >> shiftAmount());
  }

  inline int index() const {
    return int(data & ( (1 << shiftAmount()) - 1) );
  }

private:
  uint64_t data;
};


#endif /* _PR_AUXILARY_STRUCTURES_H_ */
