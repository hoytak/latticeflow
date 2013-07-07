#ifndef _NUMERICAL_H_
#define _NUMERICAL_H_

#include <cstdint>
#include <boost/multiprecision/cpp_int.hpp>

namespace latticeQBP {
  
  using namespace std;
  using boost::multiprecision::cpp_int;

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

    void add(dtype v) {
      if(uses_long_long)
        ll_sum_value += v;
      else
        inf_int_sum_value += v;

      ++n_values;
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
    void valueRoundedUp() const {
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

}; 
#endif /* _NUMERICAL_H_ */
