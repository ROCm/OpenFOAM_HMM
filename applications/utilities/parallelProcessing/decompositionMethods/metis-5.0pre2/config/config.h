
#ifndef CONFIG_H
#define CONFIG_H 1

/****************************************************************************
* A set of defines that can be modified by the user 
*****************************************************************************/
/*-------------------------------------------------------------------------- 
 Specifies the width of the elementary data type that will hold information
 about vertices and their adjacency lists. 

 Possible values:  
   32 : Use 32 bit signed integers
   64 : Use 64 bit signed integers

 A width of 64 should be specified if the number of vertices or the total
 number of edges in the graph exceed the limits of a 32 bit signed integer
 i.e., 2^31-1. 
 Proper use of 64 bit integers requires that the c99 standard datatypes 
 int32_t and int64_t are supported by the compiler.
 GCC does provides these definitions in stdint.h, but it may require some
 modifications on other architectures.
--------------------------------------------------------------------------*/
#define IDXTYPEWIDTH 32  


/*-------------------------------------------------------------------------- 
 Specifies if the __thread storage directive is available by the compiler
 to indicate thread local storage. This storage directive is available in
 most systems using gcc compiler but it may not be available in other 
 systems. 
   
 Possible values:
  0 : Not available and do not use thread local storage
  1 : It is available and the __thread modifier will be used 
--------------------------------------------------------------------------*/
#define HAVE_THREADLOCALSTORAGE 0




/****************************************************************************
* Do not change anything bellow this point
*****************************************************************************/
/* Uniform defines for various compilers */
#if defined(_MSC_VER)
  #define COMPILER_MSC
#endif
#if defined(__ICC)
  #define COMPILER_ICC
#endif
#if defined(__GNUC__)
  #define COMPILER_GCC
#endif


#if defined(COMPILER_GCC)
  #include <stdint.h>
#endif


#if defined(COMPILER_MSC)
  #include <ctrdefs.h>
  #define __thread __declspec( thread )

  typedef __int32                 int32_t;
  typedef __int64                 int64_t;
  typedef unsigned __int32        uint32_t;
  typedef unsigned __int64        uint64_t;
#endif


#if defined(UNIX)
  #include <getopt.h>
  #include <sys/time.h>
  #include <sys/resource.h>
#endif

#if defined(ENABLE_OPENMP)
  #include <omp.h>
#endif


#endif /* CONFIG_H */
