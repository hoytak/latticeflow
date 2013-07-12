#ifndef _NUMERICAL_H_
#define _NUMERICAL_H_

#include <cstdint>

#ifndef EMACS_FLYMAKE
#include <boost/multiprecision/cpp_int.hpp>
#endif

namespace latticeQBP {
  
  using namespace std;

#ifndef EMACS_FLYMAKE
  using namespace boost::multiprecision;
#endif

#ifdef EMACS_FLYMAKE

  template <typename dtype> class CompType { typedef dtype Type; };
  template <typename dtype> class DTypeMap { typedef dtype Type; };

#else

#ifndef NDEBUG
  
  template <typename dtype> class CompType {};

  template <> class CompType<int32_t> { typedef int64_t Type;          };
  template <> class CompType<int64_t> { typedef checked_int128_t Type; };

#ifdef HAVE__int128_t
  template <> class CompType<__int128_t> { typedef checked_int256_t Type; };
#endif

#else // NDEBUG defined

  template <typename dtype> class CompType {};

  template <> class CompType<int32_t> { typedef int64_t Type; };

#ifdef HAVE__int128_t
  template <> class CompType<int64_t> { typedef __int128_t Type; };
#else
  template <> class CompType<int64_t> { typedef int128_t Type; };
#endif

#ifdef HAVE__int128_t
  template <> class CompType<__int128_t> { typedef int256_t Type; };
#else
  template <> class CompType<int128_t> { typedef int256_t Type; };
#endif
  
#endif

  ////////////////////////////////////////////////////////////////////////////////
  // Now the reverse types
  
#ifndef NDEBUG
  
  template <typename comp_type> class DTypeMap {};

  template <> class DTypeMap<int64_t> { typedef int32_t Type;          };
  template <> class DTypeMap<checked_int128_t> { typedef int64_t Type; };

#ifdef HAVE__int128_t
  template <> class DTypeMap<checked_int256_t> { typedef __int128_t Type; };
#else
  template <> class DTypeMap<checked_int256_t> { typedef int128_t Type; };
#endif

#else // NDEBUG defined

  template <typename dtype> class DType {};

  template <> class DTypeMap<int64_t> { typedef int32_t Type; };

#ifdef HAVE__int128_t
  template <> class DTypeMap<__int128_t> { typedef int64_t Type; };
#else
  template <> class DTypeMap<int128_t> { typedef int64_t Type; };
#endif

#ifdef HAVE__int128_t
  template <> class DTypeMap<int256_t> { typedef __int128_t Type; };
#else
  template <> class DTypeMap<int256_t> { typedef int128_t Type; };
#endif
  
#endif

#endif

  ////////////////////////////////////////////////////////////////////////////////
  // Now casting and other low level stuff

  template <typename CompType> 
  static inline typename DTypeMap<CompType>::Type toDType(const CompType& c) {
    typedef typename DTypeMap<CompType>::Type dtype;
    
    dtype t = dtype(c);
    assert_equal(CompType(t), c);
    return t;
  }

  template <typename CompType> 
  static inline typename DTypeMap<CompType>::Type ceilAverage(CompType sum, size_t n) {
    if(sum > 0) sum += (n-1);
    return toDType(sum / n);
  }


}; 


#endif /* _NUMERICAL_H_ */

