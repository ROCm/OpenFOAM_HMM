/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * struct.h
 *
 * This file contains data structures for the mgridgen routines.
 *
 * George Irene
 */


/*************************************************************************
* This data structure holds the input graph
**************************************************************************/
struct graphdef {
  int nvtxs;            /* The # of vertices in the graph */
  idxtype *xadj;        /* Pointers to the locally stored vertices */
  idxtype *adjncy;      /* Array that stores the adjacency lists of nvtxs */
  idxtype *vwgt;        /* Vertex weights */
  realtype *vvol;       /* The volume of the vertex */
  realtype *vsurf;      /* The surface of the vertex (applicable only to boundary elements) */
  realtype *adjwgt;     /* Array that stores the weights of the adjacency lists */
  realtype *adjwgtsum;  /* The sum of the adjacent weights (surace volume) */

  idxtype *cmap;        /* Used for coarsening */

  idxtype *where;       /* Partitioning vector */
  idxtype *pwgts;       /* The weight of the partitions */
  
  int nmoves;           /* The number of moves during refinement */

  realtype *pvol;       /* The volume of the partitions */
  realtype *psurf;      /* The surface area of the partitions */

  realtype minratio;    /* The value of the objective function */

  struct graphdef *finer;
  struct graphdef *coarser;
};

typedef struct graphdef GraphType;


/*************************************************************************
* The following structure stores information used by SHC
**************************************************************************/
struct controldef {
  int dbglvl;           /* Controls the debuging output of the program */
  int CType;            /* The type of coarsening */
  int RType;            /* The type of refinement */
  int minsize;          /* The bounds on the number of elements per partition */
  int maxsize;            		
  int nparts;           /* The number of discovered partitions */
  int dim;            	/* The number of dimensions for mesh 2D or 3D */
};

typedef struct controldef CtrlType;
