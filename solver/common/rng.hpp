#ifndef _RNG_H_
#define _RNG_H_

// First get the interface out with C99 style integer objects

extern "C" {

#include <stdint.h>

  typedef uint64_t LCGState;
    
  /* Note that the state is returned by value.  This is for pure,
   * brute speed. */

  static LCGState Lcg_New(unsigned long seed) {
    // Just make sure we don't have a weak seed with integer size truncations
    return ( ( (uint64_t) seed) + 0xcb63b83e3798bbfeull);
  }
  
  static uint32_t Lcg_Next(LCGState* state) {
    (*state) = 6364136223846793005ul * (*state) + 1442695040888963407ull;
    return (uint32_t)( ((*state) >> 32) ^ (*state) );
  }
}


class FastBernoulli {
public:
  FastBernoulli(size_t seed)
    : state(Lcg_New(seed)), count_down(32), v(Lcg_Next(&state))
  {}
  
  bool operator()() {

    if( !(--count_down) ) {
      v = Lcg_Next(&state);
      count_down = 32;
    } else {
      v >>= 1;
    }
    
    return !!(v & 0x1);
  }

private:
  LCGState state;
  size_t count_down;
  uint32_t v;
};

class FastInt {
public:
  FastInt(size_t seed)
    : state(Lcg_New(seed))
  {}
  
  int operator()(int _top) {
    return Lcg_Next(&state) % _top;
  }

  int operator()(int bottom, int top) {
    return bottom + Lcg_Next(&state) % (1 + top - bottom);
  }

private:
  LCGState state;
};

class FastUInt32 {
public:
  FastUInt32(size_t seed)
    : state(Lcg_New(seed))
  {}
  
  uint32_t operator()() {
    return Lcg_Next(&state);
  }
  
private:
  LCGState state;
};


#endif /* _RNG_H_ */

