/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * sort.c
 *
 * This function contains various utility routines
 *
 * George
 */

#include "IMlib.h"

/*************************************************************************
* This function sorts an array of type KeyValueType in increasing order
**************************************************************************/
void ikeyvalsort_org(int n, IKeyValueType *nodes)
{
  qsort((void *)nodes, (size_t)n, (size_t)sizeof(IKeyValueType), IncKeyValueCmp);
}


/*************************************************************************
* This function compares 2 KeyValueType variables for sorting in inc order
**************************************************************************/
int IncKeyValueCmp(const void *v1, const void *v2)
{
  IKeyValueType *n1, *n2;

  n1 = (IKeyValueType *)v1;
  n2 = (IKeyValueType *)v2;

  return (n1->key != n2->key ? n1->key - n2->key : n1->val - n2->val);
}


/*************************************************************************
* This function sorts an array of type KeyValueType in increasing order
**************************************************************************/
void dkeyvalsort(int n, IKeyValueType *nodes)
{
  qsort((void *)nodes, (size_t)n, (size_t)sizeof(IKeyValueType), DecKeyValueCmp);
}


/*************************************************************************
* This function uses simple counting sort to return a permutation array
* corresponding to the sorted order. The keys are assumed to start from
* 0 and they are positive.  This sorting is used during matching.
**************************************************************************/
void BucketSortKeysInc(int n, idxtype max, idxtype *keys, int *tperm, int *perm)
{
  int i, ii;
  int *counts;

  counts = ismalloc(max+2, 0, "BucketSortKeysInc: counts");

  for (i=0; i<n; i++)
     counts[keys[i]]++;
  MAKECSR(i, max+1, counts);

  for (ii=0; ii<n; ii++) {
     i = tperm[ii];
     perm[counts[keys[i]]++] = i;
  }

  IMfree(&counts, LTERM);
}


/*************************************************************************
* This function compares 2 KeyValueType variables for sorting in inc order
**************************************************************************/
int DecKeyValueCmp(const void *v1, const void *v2)
{
  IKeyValueType *n1, *n2;

  n1 = (IKeyValueType *)v1;
  n2 = (IKeyValueType *)v2;

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
void RandomPermuteFine(int n, int *p, int flag)
{
  int i, u, v, tmp;

  if (flag == 1)
    for (i=0; i<n; i++)
      p[i] = i;

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
