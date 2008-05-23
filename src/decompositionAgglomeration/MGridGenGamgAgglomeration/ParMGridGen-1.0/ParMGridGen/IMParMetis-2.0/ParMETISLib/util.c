/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * util.c
 *
 * This function contains various utility routines
 *
 * George Irene
 */

#include "parmetis.h"


/*************************************************************************
* This function prints an error message and exits
**************************************************************************/
void errexit(char *f_str,...)
{
  va_list argp;
  char out1[256], out2[256];

  va_start(argp, f_str);
  vsprintf(out1, f_str, argp);
  va_end(argp);

  sprintf(out2, "Error! %s", out1);

  fprintf(stdout, out2);
  fflush(stdout);

  abort();
}


/*************************************************************************
* This function prints an error message and exits
**************************************************************************/
void myprintf(CtrlType *ctrl, char *f_str,...)
{
  va_list argp;
  char out1[256], out2[256];

  va_start(argp, f_str);
  vsprintf(out1, f_str, argp);
  va_end(argp);

  sprintf(out2, "[%2d] %s", ctrl->mype, out1);

  fprintf(stdout, out2);
  fflush(stdout);

}



/*************************************************************************
* This function prints an error message and exits
**************************************************************************/
void rprintf(CtrlType *ctrl, char *f_str,...)
{
  va_list argp;

  if (ctrl->mype == 0) {
    va_start(argp, f_str);
    vfprintf(stdout, f_str, argp);
    va_end(argp);
  }

  fflush(stdout);

  MPI_Barrier(ctrl->comm);

}


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
* The following function allocates an array of float 
**************************************************************************/
float *fmalloc(int n, char *msg)
{
  if (n == 0)
    return NULL;

  return (float *)IMmalloc(sizeof(float)*n, msg);
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
* The follwoing function allocates an array of reals
**************************************************************************/
realtype *realsmalloc(int n, realtype rval, char *msg)
{
  if (n == 0)
    return NULL;

  return realset(n, rval, (realtype *)IMmalloc(sizeof(realtype)*n, msg));
}

#endif

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
int idxamax(int n, idxtype *x)
{
  int i, max=0;

  for (i=1; i<n; i++)
    max = (x[i] > x[max] ? i : max);

  return max;
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
* This function computes a 2-norm
**************************************************************************/
float snorm2(int n, float *v)
{
  int i;
  float partial = 0;
 
  for (i = 0; i<n; i++)
    partial += v[i] * v[i];

  return sqrt(partial);
}



/*************************************************************************
* This function computes a 2-norm
**************************************************************************/
float sdot(int n, float *x, float *y)
{
  int i;
  float partial = 0;
 
  for (i = 0; i<n; i++)
    partial += x[i] * y[i];

  return partial;
}


/*************************************************************************
* This function computes a 2-norm
**************************************************************************/
void saxpy(int n, float alpha, float *x, float *y)
{
  int i;
 
  for (i=0; i<n; i++)
    y[i] += alpha*x[i];
}






/*************************************************************************
* This function sorts an array of type KeyValueType in increasing order
**************************************************************************/
void ikeyvalsort_org(int n, KeyValueType *nodes)
{
  qsort((void *)nodes, (size_t)n, (size_t)sizeof(KeyValueType), IncKeyValueCmp);
}


/*************************************************************************
* This function compares 2 KeyValueType variables for sorting in inc order
**************************************************************************/
int IncKeyValueCmp(const void *v1, const void *v2)
{
  KeyValueType *n1, *n2;

  n1 = (KeyValueType *)v1;
  n2 = (KeyValueType *)v2;

  return (n1->key != n2->key ? n1->key - n2->key : n1->val - n2->val);
}



/*************************************************************************
* This function sorts an array of type KeyValueType in increasing order
**************************************************************************/
void dkeyvalsort(int n, KeyValueType *nodes)
{
  qsort((void *)nodes, (size_t)n, (size_t)sizeof(KeyValueType), DecKeyValueCmp);
}


/*************************************************************************
* This function compares 2 KeyValueType variables for sorting in inc order
**************************************************************************/
int DecKeyValueCmp(const void *v1, const void *v2)
{
  KeyValueType *n1, *n2;

  n1 = (KeyValueType *)v1;
  n2 = (KeyValueType *)v2;

  return n2->key - n1->key;

}



/*************************************************************************
* This function does a binary search on an array for a key and returns
* the index
**************************************************************************/
int BSearch(int n, idxtype *array, int key)
{
  int a=0, b=n, c;

  while (b-a > 8) {
    c = (a+b)>>1;
    if (array[c] > key)
      b = c;
    else
      a = c;
  }

  for (c=a; c<b; c++) {
    if (array[c] == key)
      return c;
  }

  errexit("Key %d not found!\n", key);
}



/*************************************************************************
* This file randomly permutes the contents of an array.
* flag == 0, don't initialize perm
* flag == 1, set p[i] = i 
**************************************************************************/
void RandomPermute(int n, idxtype *p, int flag)
{
  int i, u, v;
  idxtype tmp;

  if (flag == 1) {
    for (i=0; i<n; i++)
      p[i] = i;
  }

  for (i=0; i<n; i++) {
    v = RandomInRange(n);
    u = RandomInRange(n);
    SWAP(p[v], p[u], tmp);
  }
}


/*************************************************************************
* This file randomly permutes the contents of an array.
* flag == 0, don't initialize perm
* flag == 1, set p[i] = i 
**************************************************************************/
void FastRandomPermute(int n, idxtype *p, int flag)
{
  int i, u, v;
  idxtype tmp;

  if (flag == 1) {
    for (i=0; i<n; i++)
      p[i] = i;
  }

  for (i=0; i<n; i+=8) {
    v = RandomInRange(n-4);
    u = RandomInRange(n-4);
    SWAP(p[v], p[u], tmp);
    SWAP(p[v+1], p[u+1], tmp);
    SWAP(p[v+2], p[u+2], tmp);
    SWAP(p[v+3], p[u+3], tmp);
  }
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
* This function returns the log2(x)
**************************************************************************/
int log2(int a)
{
  int i;

  for (i=1; a > 1; i++, a = a>>1);
  return i-1;
}
