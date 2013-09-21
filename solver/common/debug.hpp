#ifndef HK_DEBUG_H
#define HK_DEBUG_H

// #ifdef NDEBUG
// #undef NDEBUG
// #endif


#ifndef HAVE_OUTPUT
#define HAVE_OUTPUT 0
#endif

#include <iostream>
#include <cmath>

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


#if !defined(EMACS_FLYMAKE) && !defined(NDEBUG)
#warning ">>>>>>>>>>>>>>> Debug On; pass -DNDEBUG to disable <<<<<<<<<<<<<<<<<<<<<<<<<"
// #else
// #warning ">>>>>>>>>>>>>>> Debug Off <<<<<<<<<<<<<<<<<<<<<<<<<"
#endif

////////////////////////////////////////////////////////////////////////////////
// Now redefine stuff as needed.  This allows us to turn the assert
// statements on and off for various things.

#ifdef DEBUG_MODE
#undef DEBUG_MODE
#endif

#ifdef DBHERE
#undef DBHERE
#endif

#ifdef db_printval
#undef db_printval
#endif

#ifdef assert
#undef assert
#endif

#ifdef assert_equal
#undef assert_equal
#endif

#ifdef assert_notequal
#undef assert_notequal
#endif

#ifdef assert_almost_equal
#undef assert_almost_equal
#endif

#ifdef assert_geq
#undef assert_geq
#endif

#ifdef assert_gt
#undef assert_gt
#endif

#ifdef assert_leq
#undef assert_leq
#endif

#ifdef assert_lt
#undef assert_lt
#endif

#ifdef assert_close
#undef assert_close
#endif

#ifdef assert
#undef assert
#endif

#define strong_assert(x)				\
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

#define strong_assert_equal(x, y)					\
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

#define strong_assert_almost_equal(x,y)				\
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

#define strong_assert_leq(x, y)					\
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

#define strong_assert_geq(x, y)					\
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

#define strong_assert_lt(x, y)						\
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

#define strong_assert_gt(x, y)						\
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

#define strong_assert_notequal(x, y)					\
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

#define strong_assert_close(x, y, d)					\
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



#ifndef NDEBUG

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

#define assert(x)  strong_assert(x) 
#define assert_equal(x, y) strong_assert_equal(x, y)
#define assert_notequal(x, y) strong_assert_notequal(x, y)
#define assert_almost_equal(x,y)		 strong_assert_almost_equal(x,y)		
#define assert_geq(x,y) strong_assert_geq(x,y)
#define assert_gt(x,y) strong_assert_gt(x,y)
#define assert_leq(x,y) strong_assert_leq(x,y)
#define assert_lt(x,y) strong_assert_lt(x,y)
#define assert_close(x, y, d) strong_assert_close(x, y, d)
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
#define assert_close(x, y, d)

// #warning ">>>>>>>>>>>>>> Debug Off <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

#endif

