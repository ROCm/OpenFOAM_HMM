/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * memory.c
 *
 * This function contains various memory utility routines
 *
 * George Irene
 */

#include "IMlib.h"


#ifndef DMALLOC

/*************************************************************************
* The following function allocates an array of integers
**************************************************************************/
int *imalloc(int n, char *msg)
{
  if (n == 0)
    return NULL;

  return (int *)IMmalloc(sizeof(int)*n, msg);
}


/*************************************************************************
* The following function allocates an array of integers
**************************************************************************/
idxtype *idxmalloc(int n, char *msg)
{
  if (n == 0)
    return NULL;

  return (idxtype *)IMmalloc(sizeof(idxtype)*n, msg);
}


/*************************************************************************
* The following function allocates an array of double 
**************************************************************************/
double *fmalloc(int n, char *msg)
{
  if (n == 0)
    return NULL;

  return (double *)IMmalloc(sizeof(double)*n, msg);
}


/*************************************************************************
* The following function allocates an array of real
**************************************************************************/
realtype *realmalloc(int n, char *msg)
{
  if (n == 0)
    return NULL;

  return (realtype *)IMmalloc(sizeof(realtype)*n, msg);
}


/*************************************************************************
* The follwoing function allocates an array of integers
**************************************************************************/
int *ismalloc(int n, int ival, char *msg)
{
  if (n == 0)
    return NULL;

  return iset(n, ival, (int *)IMmalloc(sizeof(int)*n, msg));
}


/*************************************************************************
* The follwoing function allocates an array of integers
**************************************************************************/
idxtype *idxsmalloc(int n, idxtype ival, char *msg)
{
  if (n == 0)
    return NULL;

  return idxset(n, ival, (idxtype *)IMmalloc(sizeof(idxtype)*n, msg));
}


/*************************************************************************
* The follwoing function allocates an array of doubles
**************************************************************************/
double *fsmalloc(int n, double ival, char *msg)
{
  if (n == 0)
    return NULL;

  return fset(n, ival, (double *)IMmalloc(sizeof(double)*n, msg));
}


/*************************************************************************
* The follwoing function allocates an array of reals
**************************************************************************/
realtype *realsmalloc(int n, realtype rval, char *msg)
{
  if (n == 0)
    return NULL;

  return realset(n, rval, (realtype *)IMmalloc(sizeof(realtype)*n, msg));
}


/*************************************************************************
* This function is my wrapper around malloc
**************************************************************************/
void *IMmalloc(int nbytes, char *msg)
{
  void *ptr;

  if (nbytes == 0)
    return NULL;

  ptr = (void *)malloc(nbytes);
  if (ptr == NULL) 
    errexit("***Memory allocation failed for %s. Requested size: %d bytes", msg, nbytes);

  return ptr;
}

#endif

/*************************************************************************
* This function is my wrapper around free, allows multiple pointers    
**************************************************************************/
void IMfree(void **ptr1,...)
{
  va_list plist;
  void **ptr;

  if (*ptr1 != NULL)
    free(*ptr1);
  *ptr1 = NULL;

  va_start(plist, ptr1);

  while ((ptr = va_arg(plist, void **)) != LTERM) {
    if (*ptr != NULL)
      free(*ptr);
    *ptr = NULL;
  }

  va_end(plist);
}            
