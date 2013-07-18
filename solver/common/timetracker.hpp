#ifndef _TIMETRACKER_H_
#define _TIMETRACKER_H_

#include <ctime>
#include <sstream>
#include <cmath>
#include "debug.hpp"

class TimeTracker{
public:
  inline TimeTracker(bool start=false);
    
  inline double start();
  inline double stop();
  inline void reset();
  inline bool active() const;
  inline double elapsedSeconds() const;
  inline std::string asString() const;

  static inline std::string timeAsString(double num_seconds);
    
private:
  bool started;
  bool paused;
  clock_t start_time;
  clock_t offset;
};

inline TimeTracker::TimeTracker(bool _start)
  : started(false), paused(false), start_time(0), offset(0)
{
  if(_start) 
	start();
}

inline double TimeTracker::start()
{
  if(paused) 
    {
      start_time = clock();
      paused = false;
      return elapsedSeconds();
    }
  else
    {
      offset = 0;
      start_time = clock();
      started = true;
      return 0;
    }
}

inline void TimeTracker::reset()
{
  started = false;
  paused = false;
  offset = 0;
}

inline double TimeTracker::elapsedSeconds() const
{
  if(started)
    {
      clock_t cur_time_diff = paused ? 0 : clock() - start_time;
      double r = double(cur_time_diff + offset) / CLOCKS_PER_SEC;
      assert_geq(r, 0);
      return r;
    }
  else
    {
      return 0;
    }
}

inline double TimeTracker::stop() 
{
  if(!started)
    {
      return 0;
    }
  else if(paused)
    {
      return elapsedSeconds();
    }
  else
    {
      clock_t cur_time = clock();
      offset += cur_time - start_time;
      start_time = 0;
      paused = true;
      return elapsedSeconds();
    }
}

inline bool TimeTracker::active() const
{
  return started && !paused; 
}

inline std::string TimeTracker::asString() const
{
  return timeAsString(elapsedSeconds());
}

inline std::string TimeTracker::timeAsString(double t)
{
  std::stringstream ss;

  if(t < 0.001)
	ss << (1000000*t) << "us";
  else if(t < 1)
    ss << (1000*t) << "ms";
  else if(t < 60)
    {
      size_t nhs = size_t(floor(100*t)) % 100;
      ss << (int(floor(t))) << (nhs >= 10 ? "." : ".0") << nhs << "s";
    }
  else if(t < 3600)
	ss << int(floor(t/60)) << "m " << (int(floor(t)) % 60) << "s";
  else if(t < 86400)
	ss << int(floor(t/3600)) << "h " << ((int(floor(t)) % 3600) / 60) << "m";
  else
	ss << int(floor(t/86400)) << "d "
	   << ((int(floor(t)) % 86400) /3600) << "h "
	   << ((int(floor(t)) % 3600) / 60) << "m";

  return ss.str();
}


#endif /* _TIMETRACKER_H_ */
