/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * util.c
 *
 * This file contains utility functions
 *
 */

#include "IMlib.h"



/*************************************************************************
* This function prints an error message and exits
**************************************************************************/
void *errexit(char *f_str,...)
{
  va_list argp;

  va_start(argp, f_str);
  vfprintf(stderr, f_str, argp);
  va_end(argp);

  fprintf(stderr,"\n");
  fflush(stderr);

  exit(0);
}


/*************************************************************************
* This function returns the log2(x)
**************************************************************************/
int IMlog2(int a)
{
  int i;

  for (i=1; a > 1; i++, a = a>>1);
  return i-1;
}


/*************************************************************************
* This function returns the log2(x)
**************************************************************************/
double flog2(double a)
{
  return log(a)/log(2.0);
}


/*************************************************************************
* This function returns true if the a is a power of 2
**************************************************************************/
int ispow2(int a)
{
  for (; a%2 != 1; a = a>>1);
  return (a > 1 ? 0 : 1);
}


/*************************************************************************
* This function returns the seconds
**************************************************************************/
double seconds(void)
{
  return((double) clock()/CLOCKS_PER_SEC);
}
