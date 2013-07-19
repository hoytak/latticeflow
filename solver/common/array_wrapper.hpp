// Wrapper file to simply get the right array depending on what's
// available.

#ifndef _ARRAY_WRAPPER_H_
#define _ARRAY_WRAPPER_H_

#include <type_traits>
#include <iostream>
#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>

#include "debug.hpp"

using std::copy;
using std::copy_n;
using std::fill;

template <typename dtype> 
class Element {
public:
  Element(size_t idx, dtype v) 
    : _index(idx), _value(v)
  {}

  size_t index() const { return _index; }
  dtype value() const { return _value; }

private:
  size_t _index;
  dtype _value;
};

template <typename dtype, int _n_elements> class Array {
private:
  dtype _data[_n_elements];
public:
  Array() {}


private:
  // Set up a generic templating function method
  void inline _ResolveArray() {}

  template <typename T, typename... Args>
  void inline _ResolveArray(const T& a1, const Args&... args,
			    typename std::enable_if<std::is_convertible<T, dtype>::value>::type* = 0) {
    _data[_n_elements - sizeof...(Args) - 1] = a1;
    _ResolveArray(args...);
  }

public: 
  
  template <typename... Args>
  inline Array(const Args&... args, 
	       typename std::enable_if<sizeof...(Args) == _n_elements>::type* = 0) {
    _ResolveArray(args...);
  }

  template <typename T>
  inline Array(const T& v, typename std::enable_if<std::is_convertible<T, dtype>::value>::type* = 0) {
    fill(_data, _data + _n_elements, dtype(v));
  }
  
  template <typename FillIterator> 
  inline Array(FillIterator iter, 
	       typename std::enable_if<std::is_convertible<
		 typename std::iterator_traits<FillIterator>::value_type, dtype>
				       ::value>::type* = 0) {
    copy_n(iter, _n_elements, _data);
  }

  template <typename T>
  inline Array(const std::array<T, _n_elements>& v, 
	       typename std::enable_if<std::is_convertible<T, dtype>::value>::type* = 0)
  {
    copy(v.begin(), v.end(), _data);
  }

  template <typename T>
  inline Array(const Array<T, _n_elements>& v, 
	       typename std::enable_if<std::is_convertible<T, dtype>::value>::type* = 0)
  {
    copy(v.begin(), v.end(), _data);
  }

  inline Array(const std::array<dtype, _n_elements>& v)
  {
    copy(v.begin(), v.end(), _data);
  }

  inline Array(const Array<dtype, _n_elements>& v) 
  {
    copy(v.begin(), v.end(), _data);
  }

  inline Array(std::initializer_list<dtype> v)
  {
    copy(v.begin(), v.end(), _data);
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Setting stuff

  template <typename T>
  inline void set(std::initializer_list<T>& v, 
                  typename std::enable_if<std::is_convertible<T, dtype>::value>::type* = 0){ 
    for(size_t i = 0; i < _n_elements; ++i)
      _data[i] = dtype(v[i]);
  }

  template <typename T>
  inline const Array& operator=(const T& v) {
    fill(_data, _data + _n_elements, dtype(v));
    return *this;
  }


  ////////////////////////////////////////////////////////////////////////////////
  // Various information
  static constexpr size_t size() { return _n_elements; }

  ////////////////////////////////////////////////////////////////////////////////
  // elementwise access

  dtype& operator[](size_t idx) {
    assert_lt(idx, _n_elements);
    return _data[idx];
  }

  const dtype& operator[](size_t idx) const {
    assert_lt(idx, _n_elements);
    return _data[idx];
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Iterators

  typedef dtype* iterator;
  typedef const dtype* const_iterator;

  iterator begin() {
    return _data;
  }

  const_iterator begin() const {
    return _data;
  }

  iterator end() {
    return _data + _n_elements;
  }

  const_iterator end() const {
    return _data + _n_elements;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // 
  dtype sum() const {
    dtype e = _data[0];

    for(int i = 1; i < _n_elements; ++i)
      e += _data[1];
    
    return e;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Element components

  template <typename alt_dtype>
  inline Array operator+(const Element<alt_dtype>& x) const {
    Array out(*this);
    assert(x.index() < _n_elements);

    out[x.index()] += x.value();

    return out;
  }

  template <typename alt_dtype>
  inline Array operator-(const Element<alt_dtype>& x) const {
    Array out(*this);
    assert(x.index() < _n_elements);

    out[x.index()] -= x.value();

    return out;
  }

  template <typename alt_dtype>
  inline Array operator*(const Element<alt_dtype>& x) const {
    Array out(*this);
    assert(x.index() < _n_elements);

    out[x.index()] *= x.value();

    return out;
  }

  template <typename alt_dtype>
  inline Array operator/(const Element<alt_dtype>& x) const {
    Array out(*this);
    assert(x.index() < _n_elements);

    out[x.index()] /= x.value();

    return out;
  }


  ////////////////////////////////////////////////////////////////////////////////
  // other arrays

  template <typename alt_dtype>
  inline Array operator+(const std::array<alt_dtype, _n_elements>& x) const {
    Array out;

    for(size_t i = 0; i < _n_elements; ++i) {
      out[i] = (*this)[i] + x[i];
    }

    return out;
  }

  template <typename alt_dtype>
  Array operator-(const std::array<alt_dtype, _n_elements>& x) const {
    Array out;

    for(size_t i = 0; i < _n_elements; ++i) {
      out[i] = (*this)[i] - x[i];
    }

    return out;
  }

  template <typename alt_dtype>
  inline Array operator/(const std::array<alt_dtype, _n_elements>&  x) const {
    Array out;

    for(size_t i = 0; i < _n_elements; ++i) {
      out[i] = (*this)[i] / x[i];
    }

    return out;
  }

  template <typename alt_dtype>
  inline Array operator*(const std::array<alt_dtype, _n_elements>&  x) const {
    Array out;

    for(size_t i = 0; i < _n_elements; ++i) {
      out[i] = (*this)[i] * x[i];
    }

    return out;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Scalar components

  template <typename alt_dtype>
  inline Array operator+(const alt_dtype& x) const {
    Array out;

    for(size_t i = 0; i < _n_elements; ++i) {
      out[i] = (*this)[i] + x[i];
    }

    return out;
  }

  template <typename alt_dtype>
  Array operator-(const alt_dtype& x) const {
    Array out;

    for(size_t i = 0; i < _n_elements; ++i) {
      out[i] = (*this)[i] - x[i];
    }

    return out;
  }

  template <typename alt_dtype>
  inline Array operator/(const alt_dtype& x) const {
    Array out;

    for(size_t i = 0; i < _n_elements; ++i) {
      out[i] = (*this)[i] / x[i];
    }

    return out;
  }

  template <typename alt_dtype>
  inline Array operator*(const alt_dtype& x) const {
    Array out;

    for(size_t i = 0; i < _n_elements; ++i) {
      out[i] = (*this)[i] * x[i];
    }

    return out;
  }
}; 

template <typename T, int _n_elements> 
typename Array<T, _n_elements>::iterator begin(Array<T, _n_elements>& a) {
  return a.begin();
}

template <typename T, int _n_elements> 
typename Array<T, _n_elements>::const_iterator begin(const Array<T, _n_elements>& a) {
  return a.begin();
}

template <typename T, int _n_elements> 
typename Array<T, _n_elements>::iterator end(Array<T, _n_elements>& a) {
  return a.end();
}

template <typename T, int _n_elements> 
typename Array<T, _n_elements>::const_iterator end(const Array<T, _n_elements>& a) {
  return a.end();
}

template <typename T, int _n_elements> 
Array<T, _n_elements> abs(const Array<T, _n_elements>& a) {

  Array<T, _n_elements> ret;
  for(int i = 0; i < _n_elements; ++i)
    ret[i] = abs(a[i]);

  return ret;
}

template <typename T, int _n_elements_1, int _n_elements_2> 
Array<T, _n_elements_1 + _n_elements_2> 
cat(const Array<T, _n_elements_1>& a1, const Array<T, _n_elements_1>& a2) {

  Array<T, _n_elements_1 + _n_elements_2> ret;

  for(int i = 0; i < _n_elements_1; ++i)
    ret[i] = a1[i];

  for(int i = 0; i < _n_elements_2; ++i)
    ret[i + _n_elements_1] = a2[i];
 
  return ret;
}

template <typename T, typename Te, int _n_elements> 
Array<T, _n_elements + 1> 
cat(const Te& t, const Array<T, _n_elements>& a, 
    typename std::enable_if<std::is_convertible<Te, T>::value>::type* = 0) {

  Array<T, _n_elements + 1> ret;
  ret[0] = t;

  for(int i = 0; i < _n_elements; ++i)
    ret[i+1] = a[i];
 
  return ret;
}


template <typename T, typename Te, int _n_elements> 
Array<T, _n_elements + 1> 
cat(const Array<T, _n_elements>& a, const Te& t,
    typename std::enable_if<std::is_convertible<Te, T>::value>::type* = 0) {

  Array<T, _n_elements + 1> ret;

  for(int i = 0; i < _n_elements; ++i)
    ret[i] = a[i];
 
  ret[_n_elements] = t;

  return ret;
}


template <int start, int end, typename T, int _n_elements>
Array<T, end - start> slice(const Array<T, _n_elements>& a) {

  static_assert(end <= _n_elements, "End of slice is past end of array.");
  static_assert(start >= 0, "Start of slice is negative.");

  Array<T, end - start> ret;
  for(int i = start; i < end; ++i)
    ret[i - start] = a[i];

  return ret;
}

template <typename T1, typename T2, int _n_elements> 
double
dist2(const Array<T1, _n_elements>& a1, const Array<T2, _n_elements>& a2,
     typename std::enable_if<std::is_convertible<T1, double>::value 
			     && std::is_convertible<T2, double>::value>::type* = 0) {
  
  double d = 0;
  
  for(uint i = 0; i < _n_elements; ++i) {
    double de = (double(a1[i]) - double(a2[i]));
    d += de*de;
  }

  return d;
}

template <typename T1, typename T2, int _n_elements> 
double
dist(const Array<T1, _n_elements>& a1, const Array<T2, _n_elements>& a2,
     typename std::enable_if<std::is_convertible<T1, double>::value 
			     && std::is_convertible<T2, double>::value>::type* = 0) {
  return sqrt(dist2(a1, a2));
}

template <typename T, int _n_elements> 
T max(const Array<T, _n_elements>& a) {
  switch(_n_elements) {
  case 0: return 0;
  case 1: return a[0];
  case 2: return std::max(a[0], a[1]);
  case 3: return std::max(a[0], std::max(a[1], a[2]));
  default: return std::accumulate(a.begin() + 1, a.end(), a[0], 
				  [](const T& t1, const T& t2){return std::max(t1, t2);});
  }
}



template<class Ch,class Tr, typename T, int _n_elements>
std::basic_ostream<Ch,Tr>& operator<<(std::basic_ostream<Ch,Tr>& os, const Array<T, _n_elements>& a) {

  os << "(";
  
  for(size_t i = 0; i < _n_elements - 1; ++i)
    os << a[i] << ", ";
    
  os << a[_n_elements - 1] << ")";

  return os;
}

template <typename Arg1, typename... Args>
inline Array<Arg1, sizeof...(Args) + 1> Ar(const Arg1& a1, const Args&... args) {
  return Array<Arg1, sizeof...(Args) + 1>({a1, args...});
}



#endif /* _ARRAY_DIRECTOR_H_ */
