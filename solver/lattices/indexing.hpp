#ifndef _INDEXING_H_
#define _INDEXING_H_

#include "../common.hpp"

class Neighbor {
public:
  Neighbor(int __dim, int __direction) 
    : _dim(__dim), _direction( direction_type(__direction))
  {
    assert(__direction == -1 || __direction == 1);
  }

  bool operator==(const Neighbor& n) const {
    return (_direction == n._direction) && (_dim == n._dim);
  }

  bool operator!=(const Neighbor& n) const {
    return (_direction != n._direction) || (_dim != n._dim);
 }

  int dim() const { return _dim; }
  int direction() const {return int(_direction); }
  int reversed() const {return -int(_direction); }
					
private:
  typedef enum {Back = -1, Forward = 1} direction_type;
  int _dim;
  direction_type _direction;
};

template <int nd>
class Indexer {
public:
  
  typedef Array<size_t, nd> index_vect;

  Indexer(const index_vect& __dimensions)
    : _dimensions(__dimensions)
    , _stride(dim_vector_factors(_dimensions))
  {
    assert_equal(prod(_dimensions), size());
  }

  Indexer(Indexer&&) = default;
  
  long neighborStep(int dim, int step) const {
    return step * _stride[dim];
  }

  template <typename step_type>
  long neighborStep(const Array<step_type, nd>& diff) const {
    long jump = 0;

    for(size_t dim = 0; dim < nd; ++dim) {
      jump += diff[dim] * _stride[dim];
    }

    return jump;
  }

  size_t neighborIndex(size_t index, int dim, int step) const {
    return index + neighborStep(dim, step);
  }

  size_t neighborIndex(const index_vect& coords, int dim, int step) const {
    return getIndex(coords) + neighborStep(dim, step);
  }

  size_t neighborIndex(size_t index, const Neighbor& n) const {
    return index + neighborStep(n.dim(), n.direction());
  }

  size_t neighborIndex(const index_vect& coords, const Neighbor& n) const {
    return getIndex(coords) + neighborStep(n.dim(), n.direction());
  }
  
  size_t getIndex(const index_vect& coords) const {
    size_t r_idx = 0;

    for(size_t i = 0; i < nd; ++i)
      r_idx += coords[i] * _stride[i];

    return r_idx;
  }

  size_t operator[](const index_vect& iv) const {
    return getIndex(iv);
  }

  index_vect getCoords(size_t index) const {

    index_vect coords;
    size_t _index = index;

    for(size_t i = nd-1; i != 0; --i) {
      coords[i] = _index % _dimensions[i];
      _index /= _dimensions[i];
    }
    coords[0] = _index;

    assert_equal(getIndex(coords), index);
    
    return coords;
  }

  index_vect operator[](size_t index) const {
    return getCoords(index);
  }

  const index_vect& dimensions() const {
    return _dimensions;
  }

  size_t stride(int dim) const {
    return _stride[dim];
  }

  size_t size() const {
    return _stride[0]*_dimensions[0];
  }

  size_t size(size_t dim) const {
    return _dimensions[dim];
  }

  unsigned int indexDimension(int idx) const {
    return idx % nd;
  }

  unsigned int indexDirection(int idx) const {
    return (idx < nd) ? -1 : 1;
  }

private:
  index_vect _dimensions;
  index_vect _stride;
};

////////////////////////////////////////////////////////////////////////////////
//


template <int nd>
class IndexIterator {
public:
  typedef typename Indexer<nd>::index_vect index_vect;
  typedef Indexer<nd> Indexer_type;

  IndexIterator(const index_vect& __dimensions = index_vect(0) ) 
    : _index(0)
    , _coords(0)
    , _dimensions(__dimensions)
  {}

  // IndexIterator(const Indexer_type& indexer)
  //   : _index(0)
  //   , end_index(prod(indexer.dimensions()))
  //   , _coords(0)
  //   , _dimensions(indexer.dimensions())
  // {}

  IndexIterator& operator++() {
    ++_index;

    for(size_t i = nd; i != 0;) {
      --i;
      ++_coords[i];

      if(i != 0 &&_coords[i] == _dimensions[i]) {
	_coords[i] = 0;
	continue;
      } else {
	break;
      }
    }

    assert_equal(Indexer<nd>(_dimensions).getIndex(_coords), _index);

    return *this;
  }

  void jump(long increment) {
    _index += increment;
    _coords = Indexer<nd>(_dimensions).getCoords(_index);

    // cout << "_coords = " << _coords << "; index = " << _index << endl;

    assert_equal(Indexer<nd>(_dimensions).getIndex(_coords), _index);
  }

  bool done() const {
    return _coords[0] >= _dimensions[0];
  }

  const index_vect& coords() const {
    return _coords;
  }

  size_t index() const {
    return _index;
  }

  size_t operator[](size_t i) const {
    return _coords[i];
  }

private:
  size_t _index;
  index_vect _coords;
  index_vect _dimensions;
  bool _done;
};

template <int nd>
class BoundedIndexIterator {
public:
  typedef typename Indexer<nd>::index_vect index_vect;

  BoundedIndexIterator(const index_vect& __dimensions = index_vect(0), 
		       const index_vect& __bounds = index_vect(0) ) 
    : _index(0)
    , indexer(__dimensions)
    , _bounds(__bounds)
    , _coords(0)
    , _bounded_index(0)
  {}

  inline size_t increment() {
    ++_bounded_index;
    ++_index;

    bool reset_index = false;

    for(size_t i = nd; i != 0;) {
      --i;
      ++_coords[i];

      if(i != 0 && _coords[i] == _bounds[i]) {
	_coords[i] = 0;
	reset_index = true;
	continue;
      } else {
	break;
      }
    }

    size_t ret = 1;
    
    if(reset_index) {
      size_t old_index = _index - 1;
      _index = indexer.getIndex(_coords);
      ret = _index - old_index;
    }

    // Skip 0, as that's the done indicator
    for(size_t i = 1; i < nd; ++i)
      assert_lt(_coords[i], _bounds[i]);

    assert_equal(_index, indexer.getIndex(_coords));

    return ret;
  }

  BoundedIndexIterator& operator++() {
    increment();
    return *this;
  }

  bool done() const {
    return _coords[0] >= _bounds[0];
  }

  const index_vect& coords() const {
    return _coords;
  }

  size_t index() const {
    return _index;
  }

  size_t boundedIndex() const {
    return _bounded_index;
  }

  size_t operator[](size_t i) const {
    return _coords[i];
  }

private:
  size_t _index;
  Indexer<nd> indexer;
  index_vect _bounds;
  index_vect _coords;
  bool _done;
  size_t _bounded_index;
};


#endif /* _INDEXING_H_ */
