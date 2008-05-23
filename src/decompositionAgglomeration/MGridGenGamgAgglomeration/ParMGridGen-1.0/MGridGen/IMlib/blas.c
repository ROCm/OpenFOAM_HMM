/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * blas.c
 *
 * This function contains various utility routines
 *
 * George Irene
 */

#include "IMlib.h"


/*************************************************************************
* These functions set the values of a vector
**************************************************************************/
int *iset(int n, int val, int *x)
{
  int i;

  for (i=0; i<n; i++)
    x[i] = val;

  return x;
}


/*************************************************************************
* These functions set the values of a vector
**************************************************************************/
idxtype *idxset(int n, idxtype val, idxtype *x)
{
  int i;

  for (i=0; i<n; i++)
    x[i] = val;

  return x;
}


/*************************************************************************
* These functions set the values of a vector
**************************************************************************/
double *fset(int n, double val, double *x)
{
  int i;

  for (i=0; i<n; i++)
    x[i] = val;

  return x;
}


/*************************************************************************
* These functions set the values of a vector
**************************************************************************/
realtype *realset(int n, realtype val, realtype *x)
{
  int i;

  for (i=0; i<n; i++)
    x[i] = val;

  return x;
}


/*************************************************************************
* These functions return the index of the maximum element in a vector
**************************************************************************/
int iamax(int n, int *x)
{
  int i, max=0;

  for (i=1; i<n; i++)
    max = (x[i] > x[max] ? i : max);

  return max;
}


/*************************************************************************
* These functions return the index of the maximum element in a vector
**************************************************************************/
int idxamax(int n, idxtype *x)
{
  int i, max=0;

  for (i=1; i<n; i++)
    max = (x[i] > x[max] ? i : max);

  return max;
}


/*************************************************************************
* These functions return the index of the maximum element in a vector
**************************************************************************/
int famax(int n, double *x)
{
  int i, max=0;

  for (i=1; i<n; i++)
    max = (x[i] > x[max] ? i : max);

  return max;
}


/*************************************************************************
* These functions return the index of the minimum element in a vector
**************************************************************************/
int iamin(int n, int *x)
{
  int i, min=0;

  for (i=1; i<n; i++)
    min = (x[i] < x[min] ? i : min);

  return min;
}


/*************************************************************************
* These functions return the index of the minimum element in a vector
**************************************************************************/
int idxamin(int n, idxtype *x)
{
  int i, min=0;

  for (i=1; i<n; i++)
    min = (x[i] < x[min] ? i : min);

  return min;
}


/*************************************************************************
* These functions return the index of the minimum element in a vector
**************************************************************************/
int famin(int n, double *x)
{
  int i, min=0;

  for (i=1; i<n; i++)
    min = (x[i] < x[min] ? i : min);

  return min;
}


/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
int charsum(int n, char *x)
{
  int i, sum = 0;

  for (i=0; i<n; i++)
    sum += x[i];

  return sum;
}


/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
int isum(int n, int *x)
{
  int i, sum = 0;

  for (i=0; i<n; i++)
    sum += x[i];

  return sum;
}


/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
int idxsum(int n, idxtype *x)
{
  int i, sum = 0;

  for (i=0; i<n; i++)
    sum += x[i];

  return sum;
}


/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
double ssum(int n, double *x)
{
  int i;
  double sum = 0.0;

  for (i=0; i<n; i++)
    sum += x[i];

  return sum;
}


/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
double ssum_strd(int n, double *x, int incx)
{
  int i;
  double sum = 0.0;

  for (i=0; i<n; i++, x+=incx)
    sum += *x;

  return sum;
}


/*************************************************************************
* This function scales the entries in an array
**************************************************************************/
void sscale(int n, double alpha, double *x)
{
  int i;

  for (i=0; i<n; i++)
    x[i] *= alpha;
}


/*************************************************************************
* This function computes a 2-norm
**************************************************************************/
double snorm2(int n, double *v)
{
  int i;
  double partial = 0;
 
  for (i = 0; i<n; i++)
    partial += v[i] * v[i];

  return sqrt(partial);
}


/*************************************************************************
* This function computes a dot product
**************************************************************************/
double sdot(int n, double *x, double *y)
{
  int i;
  double partial = 0;
 
  for (i = 0; i<n; i++)
    partial += x[i] * y[i];

  return partial;
}


/*************************************************************************
* This function computes a saxpy operation
**************************************************************************/
void saxpy(int n, double alpha, double *x, int incx, double *y, int incy)
{
  int i;
 
  for (i=0; i<n; i++, x+=incx, y+=incy) 
    *y += alpha*(*x);
}
