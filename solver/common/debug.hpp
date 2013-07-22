#ifndef HK_DEBUG_H
#define HK_DEBUG_H

#include <iostream>
#include <cmath>

#ifdef assert
#undef assert
#endif

#ifndef NDEBUG

#ifndef EMACS_FLYMAKE
#warning ">>>>>>>>>>>>>>> Debug On; pass -DNDEBUG to disable <<<<<<<<<<<<<<<<<<<<<<<<<"
#endif

#define DEBUG_MODE true

#define DBHERE					\
  do{						\
    std::cout << __FILE__ << ":"		\
	      << __FUNCTION__ << ":"		\
	      << __LINE__ << ": "		\
	      << "HERE" << std::endl;		\
  }while(0)
    

#define db_printval(x)					\
  do{							\
    std::cout << __FILE__ << ":"			\
	      << __FUNCTION__ << ":"			\
	      << __LINE__ << ": "			\
	      << #x << " = " << (x) << std::endl;	\
  }while(0)

#define assert(x)				\
  do{						\
    if(!(x))					\
      {						\
	std::cout << "ASSERTION FAILED: "	\
		  << __FILE__ << ":"		\
		  << __FUNCTION__ << ":"	\
		  << __LINE__ << ": "		\
		  << #x << std::endl;		\
	abort();				\
      }						\
  }while(0)

#define assert_equal(x, y)					\
  do{								\
    if( (x) != (y) )						\
      {								\
	std::cout << "ASSERT == FAILED: "			\
		  << __FILE__ << ":"				\
		  << __FUNCTION__ << ":"			\
		  << __LINE__ << ": "				\
		  << "(" << #x << ") " << (x)			\
		  << " != "					\
		  << (y) << " (" << #y << ") " << std::endl;	\
	abort();						\
      }								\
  }while(0)

#define assert_almost_equal(x,y)				\
  do{								\
    double xv = (x);						\
    double yv = (y);						\
    if( abs(yv - xv) > max(abs(xv), abs(yv))*1e-4 )		\
      {								\
	std::cout << "ASSERT == FAILED: "			\
		  << __FILE__ << ":"				\
		  << __FUNCTION__ << ":"			\
		  << __LINE__ << ": "				\
		  << "(" << #x << ") " << (x)			\
		  << " != "					\
		  << (y) << " (" << #y << ") " << std::endl;	\
	abort();						\
      }								\
  }while(0)

#define assert_leq(x, y)					\
  do{								\
    if(!( (x) <= (y) ))						\
      {								\
	std::cout << "ASSERT <= FAILED: "			\
		  << __FILE__ << ":"				\
		  << __FUNCTION__ << ":"			\
		  << __LINE__ << ": "				\
		  << "(" << #x << ") " << (x)			\
		  << " > "					\
		  << (y) << " (" << #y << ") " << std::endl;	\
	abort();						\
      }								\
  }while(0)

#define assert_geq(x, y)					\
  do{								\
    if(!( (x) >= (y) ))						\
      {								\
	std::cout << "ASSERT >= FAILED: "			\
		  << __FILE__ << ":"				\
		  << __FUNCTION__ << ":"			\
		  << __LINE__ << ": "				\
		  << "(" << #x << ") " << (x)			\
		  << " < "					\
		  << (y) << " (" << #y << ") " << std::endl;	\
	abort();						\
      }								\
  }while(0)

#define assert_lt(x, y)						\
  do{								\
    if(!( (x) < (y) ))						\
      {								\
	std::cout << "ASSERT < FAILED: "			\
		  << __FILE__ << ":"				\
		  << __FUNCTION__ << ":"			\
		  << __LINE__ << ": "				\
		  << "(" << #x << ") " << (x)			\
		  << " >= "					\
		  << (y) << " (" << #y << ") " << std::endl;	\
	abort();						\
      }								\
  }while(0)

#define assert_gt(x, y)						\
  do{								\
    if(!( (x) > (y) ))						\
      {								\
	std::cout << "ASSERT > FAILED: "			\
		  << __FILE__ << ":"				\
		  << __FUNCTION__ << ":"			\
		  << __LINE__ << ": "				\
		  << "(" << #x << ") " << (x)			\
		  << " >= "					\
		  << (y) << " (" << #y << ") " << std::endl;	\
	abort();						\
      }								\
  }while(0)

#define assert_notequal(x, y)					\
  do{								\
    if( (x) == (y) )						\
      {								\
	std::cout << "ASSERT != FAILED: "			\
		  << __FILE__ << ":"				\
		  << __FUNCTION__ << ":"			\
		  << __LINE__ << ": "				\
		  << "(" << #x << ") " << (x)			\
		  << " == "					\
		  << (y) << " (" << #y << ") " << std::endl;	\
	abort();						\
      }								\
  }while(0)

#define assert_close(x, y, d)					\
  do{								\
    if( abs((x) - (y)) > (d))                                   \
      {                                                                 \
	std::cout << "ASSERT =~= FAILED (d = " << (d) << "): "          \
		  << __FILE__ << ":"				\
		  << __FUNCTION__ << ":"			\
		  << __LINE__ << ": "				\
		  << "(" << #x << ") " << (x)			\
		  << " =~/~= "					\
		  << (y) << " (" << #y << ") " << std::endl;	\
	abort();						\
      }								\
  }while(0)

#else

#define DEBUG_MODE false

#define DBHERE
#define db_printval(x)
#define assert(x) 
#define assert_equal(x, y)
#define assert_notequal(x, y)
#define assert_almost_equal(x,y)		
#define assert_geq(x,y)
#define assert_gt(x,y)
#define assert_leq(x,y)
#define assert_lt(x,y)

// #warning ">>>>>>>>>>>>>> Debug Off <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

#endif

// A few overloaded things to help with debugging.

#include <iostream>
#include <vector>

template <typename T> 
std::basic_ostream<char>& operator<<(std::basic_ostream<char>& out, const std::vector<T>& v) {
  out << '(';
  for(const T& t : v)
    out << t << '|';
  out << ')' << std::endl;
  return out;
}

#endif
