#ifndef _NUMERICAL_H_
#define _NUMERICAL_H_

#include <cstdint>
#include <boost/multiprecision/cpp_int.hpp>

namespace latticeQBP {
  
  using namespace std;
  using namespace boost::multiprecision;

  ////////////////////////////////////////////////////////////////
  // This could be optimized more.  Now 
  template <typename dtype> class StableAverage
  {
  public:
    StableAverage()
      : n_values(0)
      , inf_int_sum_value(0)
      , ll_sum_value(0)
    {}

    static constexpr bool uses_long_long = (sizeof(long long) > sizeof(dtype)); 

    inline void add(dtype v, size_t count = 1) {
      if(uses_long_long)
        ll_sum_value += v;
      else
        inf_int_sum_value += v;

      n_values += count;
    }

    void value() const {
      if(uses_long_long)
        return dtype(ll_sum_value / n_values);
      else
        return dtype(inf_int_sum_value / n_values);
    }

  private:
    template <typename int_type> 
    inline dtype _valueRoundedUp(int_type& v) {
      // Try to avoid any copy operations of the inf int type

      if(v > 0)
        v += (n_values - 1);

      dtype average = value();

      if(v > 0)
        v -= (n_values - 1);

      return average;
    }

  public:
    inline void valueRoundedUp() const {
      if(uses_long_long)
        return _valueRoundedUp(ll_sum_value);
      else
        return _valueRoundedUp(inf_int_sum_value);
    }

  private:
    cpp_int inf_int_sum_value;
    long long ll_sum_value;
    size_t n_values;
  };

#ifndef NDEBUG
  
  template <typename dtype> class CompType {};

  template <> class CompType<int32_t> {
    typedef int64_t Type;
  };

  template <> class CompType<int64_t> {
    typedef checked_int128_t Type;
  };

#ifdef HAVE__int128_t
  template <> class CompType<__int128_t> {
    typedef checked_int256_t Type;
  };
#endif

#else

  template <typename dtype> class CompType {};

  template <> class CompType<int32_t> {
    typedef int64_t Type;
  };

#ifdef HAVE__int128_t
  template <> class CompType<int64_t> {
    typedef __int128_t Type;
  };
#else
  template <> class CompType<int64_t> {
    typedef int128_t Type;
  };
#endif

#ifdef HAVE__int128_t
  template <> class CompType<__int128_t> {
    typedef int256_t Type;
  };
#endif
  
#endif

}; 

#endif /* _NUMERICAL_H_ */

