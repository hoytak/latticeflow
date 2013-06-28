#ifndef _TRACKING_H_
#define _TRACKING_H_

#include <google/dense_hash_map>
#include <stdint.h>

const static uint32_t* _hash_k1_lookup = NULL;
const static uint32_t* _hash_k2_lookup = NULL;
static size_t _hash_lookup_size = 0;

template <typename Value> class StateTracker {

private:
  ////////////////////////////////////////
  // Stuff for the state tracking / hashing

  struct StateKey {
    StateKey(uint32_t _k1, uint32_t _k2) 
      : k1(_k1), k2(_k2) 
    {}

    StateKey()
      : k1(uint32_t(-1)), k2(uint32_t(-1))
    {}

    uint32_t k1, k2;
  };

  size_t n;
  const StateKey empty_key;
  StateKey current_key;

  struct StateKeyHash {
    uint32_t operator()(const StateKey& k) const {
      return k.k1;
    }
  };

  struct StateKeyEqual {
    bool operator()(const StateKey& k1, const StateKey& k2) const {
      return k1.k2 == k2.k2;
    }
  };
  
  typedef google::dense_hash_map<StateKey, Value, StateKeyHash, StateKeyEqual> hm_type;

  hm_type hm;

public:
  StateTracker(size_t _n, size_t initial_size) 
    : n(_n) 
    , empty_key()
    , current_key(0,0)
    , hm(initial_size)
  {
    hm.set_empty_key(empty_key);
    refreshHashKeyLookups();
  }

  void resetState() {
    current_key = StateKey(0,0);
  }

  void flip(size_t i) {
    current_key.k1 ^= _hash_k1_lookup[i];
    current_key.k2 ^= _hash_k2_lookup[i];
  }

  template <class ArrayType>
  void setToState(const ArrayType& state) {
    resetState();
    for(size_t i = 0; i < n; ++i) {
      if(state[i])
	flip(i);
    }
  }

  void registerCurrentState(const Value& v) {
    hm.insert(make_pair(current_key, v));
  }

  Value& value() {
    return hm[current_key];
  }
  
  const Value* valueOfCurrentState() const {
    typename hm_type::const_iterator iterator = hm.find(current_key);
    if(iterator == hm.end())
      return NULL;
    else
      return &(iterator->second);
  }

  size_t currentHash() const {
    return current_key.k1;
  }

private:

  void refreshHashKeyLookups() {
    
    if(_hash_lookup_size < n) {

#ifdef USE_OPENMP 
#pragma omp critical
#endif
      {
	if(_hash_lookup_size < n) {
	  uint32_t *new_k1_values = new uint32_t[n];
	  uint32_t *new_k2_values = new uint32_t[n];
	  
	  FastUInt32 k1_gen(0);
	  FastUInt32 k2_gen(1);

	  for(size_t i = 0; i < n; ++i) {
	    new_k1_values[i] = k1_gen();
	    new_k2_values[i] = k2_gen();
	  }

	  // Done in stages for thread safety for reading these
	  // values.
	  const uint32_t *old_k1_values = _hash_k1_lookup;
	  const uint32_t *old_k2_values = _hash_k2_lookup;
	  _hash_k1_lookup = new_k1_values;
	  _hash_k2_lookup = new_k2_values;
	  _hash_lookup_size = n;

	  if(old_k1_values == NULL)
	    delete old_k1_values;
	  if(old_k2_values == NULL)
	    delete old_k2_values;
	}
      }
    }
  }
};

#endif /* _TRACKING_H_ */


