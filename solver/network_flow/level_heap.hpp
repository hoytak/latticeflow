#ifndef _PR_EXPANSION_LEVEL_HEAP_H_
#define _PR_EXPANSION_LEVEL_HEAP_H_

#include <vector>
#include "../common.hpp"

namespace latticeQBP {

  using std::vector;

  template <class T> class PRLevelHeap {
  public:

    PRLevelHeap(size_t n_reserve_levels) 
      : levels(n_reserve_levels)
      , cur_level(~size_t(0))
      , num_levels(0)
      , count(0)
    {}

    inline size_t size() const {
      return count;
    }

    inline bool empty() const {
      return count == 0;
    }

    void reset() {
      count = 0;
      num_levels = 0;
      cur_level = ~size_t(0);

      for(auto& v : levels)
        v.clear();
    }
  
    inline void push(const size_t level, const T& t) {
    
      if(cur_level != ~size_t(0))
        assert_geq(level, cur_level);

      if(unlikely(level >= levels.size())) {
        levels.resize(2*level + 1);
      }

      bool reset_iter = (cur_level == level);

      size_t index = (reset_iter
                      ? top_iterator - levels[cur_level].begin()
                      : 0);

      levels[level].push_back(t);
    
      if(reset_iter)
        top_iterator = levels[cur_level].begin() + index;

      ++count;

      if(num_levels <= level) {
        num_levels = level + 1;
      }
    }

    inline void pushCurrentLevel(const T& t) {

      assert_lt(cur_level, levels.size());

      size_t index = top_iterator - levels[cur_level].begin();

      levels[cur_level].push_back(t);
      ++count;
    
      top_iterator = levels[cur_level].begin() + index;
    }

    // Must be called before top or pop
    inline void initIteration() {
    
      cur_level = 0;

      while(levels[cur_level].empty()) {
        assert(cur_level != num_levels);
        ++cur_level;
      }
    
      top_iterator = levels[cur_level].begin();
    }

    inline const T& current() const {
      assert(top_iterator != levels[cur_level].end());
      assert(!empty());

      return *top_iterator;
    }

    inline size_t currentLevel() const {
      return cur_level;
    }

    inline void pop() {
      assert(!empty());

      ++top_iterator;
      --count;

      if(unlikely(top_iterator == levels[cur_level].end())) {
        levels[cur_level].clear();

        if(count == 0) {
          return;
        }

        do{
          ++cur_level;
          assert_lt(cur_level, num_levels);
        } while(levels[cur_level].empty());
      
        top_iterator = levels[cur_level].begin();
      }
    }

  private:
    typedef typename vector<T>::const_iterator iterator;

    vector<vector<T> > levels;

    iterator top_iterator;
    size_t cur_level;
    size_t num_levels;
    size_t count;
  };

}; 
#endif /* _LEVEL_DEQUE_H_ */
