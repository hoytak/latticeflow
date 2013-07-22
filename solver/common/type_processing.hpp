#ifndef _TYPE_PROCESSING_H_
#define _TYPE_PROCESSING_H_

#include <iostream>
#include <cstdint>
#include <boost/multiprecision/cpp_int.hpp>

namespace latticeQBP {

  using namespace std;

  typedef __int128 int128_type;

#ifndef NDEBUG
  typedef boost::multiprecision::checked_int256_t int256_type;
#else
  typedef boost::multiprecision::int256_t int256_type;
#endif

  template <typename dtype> struct CompType {};
  template <> struct CompType<int32_t> { typedef int64_t Type; };
  template <> struct CompType<int64_t> { typedef int128_type Type; };
  template <> struct CompType<int128_type> { typedef int256_type Type; };

  ////////////////////////////////////////////////////////////////////////////////
  // Now the reverse types
  template <typename comp_type> struct DTypeMap {};
  template <> struct DTypeMap<int64_t> { typedef int32_t Type; };
  template <> struct DTypeMap<int128_type> { typedef int64_t Type; };
  template <> struct DTypeMap<int256_type> { typedef int128_type Type; };

  ////////////////////////////////////////////////////////////////////////////////
  // Now casting and other low level stuff

  template <typename CompType> 
  static inline typename DTypeMap<CompType>::Type toDType(const CompType& c) {
    typedef typename DTypeMap<CompType>::Type dtype;
    
    dtype t = static_cast<dtype>(c);
    assert_equal(static_cast<CompType>(t), c);
    return t;
  }

  template <typename CompType> 
  static inline typename DTypeMap<CompType>::Type ceilAverage(CompType sum, long n) {
    if(sum > 0) sum += (n-1);
    sum /= n;
    return toDType(sum);
  }

  template <typename CompType> 
  static inline typename DTypeMap<CompType>::Type floorAverage(CompType sum, long n) {
    if(sum < 0) sum -= (n-1);
    sum /= n;
    return toDType(sum);
  }
}; 

std::basic_ostream<char>& operator<<(std::basic_ostream<char>& out, const latticeQBP::int128_type& t) {
  out << int64_t(t / (latticeQBP::int128_type(1) << 64)) << ':' << int64_t(t & (~int64_t(0)));
  return out;
}

#endif /* _NUMERICAL_H_ */

