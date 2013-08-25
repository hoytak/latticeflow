#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include <algorithm>
#include <vector>
#include "array_wrapper.hpp"

template <int nd>
Array<size_t, nd> dim_vector_factors(const Array<size_t, nd>& a) {

  Array<size_t, nd> dim_factor;

  // Make the dim factors
  dim_factor[nd - 1] = 1;

  for(size_t i = nd - 1; i != 0; --i)
    dim_factor[i - 1] = dim_factor[i] * a[i];
  
  return dim_factor;
}

template <typename dtype, typename T, int nd>
Array<dtype, nd + 1> concat(const Array<dtype, nd>& a, const T& v,
                            typename std::enable_if<std::is_convertible<T, dtype>::value>::type* = 0) {

  Array<dtype, nd+1> ret;

  for(size_t i = 0; i < nd; ++i)
    ret[i] = a[i];

  ret[nd] = static_cast<dtype>(v);

  return ret;
}

template <typename dtype, typename T, int nd>
Array<dtype, nd + 1> concat(T v, const Array<dtype, nd>& a,
                            typename std::enable_if<std::is_convertible<T, dtype>::value>::type* = 0) {
  Array<dtype, nd+1> ret;

  ret[0] = static_cast<dtype>(v);

  for(size_t i = 0; i < nd; ++i)
    ret[i+1] = a[i];

  return ret;
}

template <typename dtype, int nd1, typename T, int nd2>
Array<dtype, nd1 + nd2> concat(const Array<dtype, nd1>& a1, const Array<T, nd2>& a2, 
                               typename std::enable_if<std::is_convertible<T, dtype>::value>::type* = 0) {

  Array<dtype, nd1 + nd2> ret;

  for(int i = 0; i < nd1; ++i)
    ret[i] = a1[i];

  for(int i = 0; i < nd2; ++i)
    ret[i+nd1] = static_cast<dtype>(a2[i]);

  return ret;
}

template <typename dtype, int nd>
Array<dtype, nd> add_each(const Array<dtype, nd>& a, dtype v) {

  Array<dtype, nd> ret;

  for(size_t i = 0; i < nd; ++i)
    ret[i] = a[i] + v;

  return ret;
}

template <typename dtype, int d>
Array<dtype, d> add_dim(const Array<dtype, d>& a, size_t idx, dtype v) {

  Array<dtype, d> ret;

  for(size_t i = 0; i < d; ++i)
    ret[i] = a[i];

  ret[idx] += v;

  return ret;
}

template <typename dtype, int d>
dtype prod(const Array<dtype, d>& a) {
  dtype x = 1;
  for(size_t i = 0; i < d; ++i)
    x *= a[i];
  return x;
}

template <typename dtype>
dtype prod(const std::vector<dtype>& a) {
  dtype x = 1;
  for(size_t i = 0; i < a.size(); ++i)
    x *= a[i];
  return x;
}


template <typename dtype, int d>
dtype sum(const Array<dtype, d>& a) {
  dtype x = 0;
  for(size_t i = 0; i < d; ++i)
    x += a[i];
  return x;
}

template <typename dtype>
Array<dtype, 1> make_array_1(const dtype& v0) {
  Array<dtype, 1> a;
  a[0] = v0;
  return a;
}

template <typename dtype>
Array<dtype, 2> make_array_2(const dtype& v0, const dtype& v1) {
  Array<dtype, 2> a;
  a[0] = v0;
  a[1] = v1;
  return a;
}

template <typename dtype>
Array<dtype, 3> make_array_3(const dtype& v0, const dtype& v1, const dtype& v2) {
  Array<dtype, 3> a;
  a[0] = v0;
  a[1] = v1;
  a[2] = v2;
  return a;
}

#endif /* _UTILITIES_H_ */
