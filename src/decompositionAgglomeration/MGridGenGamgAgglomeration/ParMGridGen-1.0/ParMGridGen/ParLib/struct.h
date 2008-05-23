/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * struct.h
 *
 * This file contains data structures.
 *
 * George Irene
 */

/* Indexes are as long as integers for now */
#ifdef IDXTYPE_INT
#define IDX_DATATYPE	MPI_INT
#else
#define IDX_DATATYPE	MPI_SHORT
#endif

/* floats for now */
#ifdef TYPE_REAL
#define REAL_DATATYPE	MPI_FLOAT
#else
#define REAL_DATATYPE	MPI_DOUBLE
#endif


/*************************************************************************
* The following data structure stores key-value pair
**************************************************************************/
struct KeyValueType {
  idxtype key;
  idxtype val;
};

typedef struct KeyValueType KeyValueType;


/*************************************************************************
* The following data structure is used to store the buckets for the 
* refinment algorithms
**************************************************************************/
struct PQueueType {
  int nnodes;
  int maxnnodes;
  idxtype *perm, *iperm, *values;  
  /* iperm[i] stores where the ith entry is located
     perm[i] stores the entry that is located in the ith position */
};

typedef struct PQueueType PQueueType;


/*************************************************************************
* The following data structure stores an edge
**************************************************************************/
struct edgedef {
  realtype ewgt;
  idxtype edge;
};
typedef struct edgedef EdgeType;


/*************************************************************************
* This data structure holds various working space data
**************************************************************************/
struct workspacedef {
  idxtype *core;			/* Where pairs, indices, and degrees are coming from */
  int maxcore;

  int nlarge;				/* The size of 'Large' */

  KeyValueType *pairs;			/* Large pair array used during setup */
  idxtype *indices;			/* Large array of indxtype used for various purposes */

  /* Auxiliary parameters */
  idxtype *pv1, *pv2, *pv4;	/* Vectors of npes+1 size used in various places */
  KeyValueType *pepairs1, *pepairs2;

  EdgeType *degrees;
};

typedef struct workspacedef MGridWorkSpaceType;


/*************************************************************************
* The following data structure holds information on degrees for k-way
* partition
**************************************************************************/
struct rinfodef {
 realtype id, ed;            /* ID/ED of edges */
 int ndegrees;          /* The number of different ext-degrees */
 EdgeType *degrees;     /* List of edges */
};

typedef struct rinfodef RInfoType;


/*************************************************************************
* The following data structure holds information on degrees for k-way
* partition
**************************************************************************/
struct nrinfodef {
 int edegrees[2];  
};

typedef struct nrinfodef NRInfoType;


/*************************************************************************
* This data structure holds the input graph
**************************************************************************/
struct graphdef {
  int gnvtxs, nvtxs, nedges;
  idxtype *xadj;	/* Pointers to the locally stored vertices */
  idxtype *vwgt;	/* Vertex weights */
  realtype *vvol;	/* The volume of the vertex */
  realtype *vsurf;	/* The surface of the vertex (applicable only to boundary elements) */
  idxtype *vsize;	/* Vertex size */
  idxtype *adjncy;	/* Array that stores the adjacency lists of nvtxs */
  realtype *adjwgt;	/* Array that stores the weights of the adjacency lists */
  realtype *adjwgtsum;	/* The sum of the adjacent weights (surace volume) */
  idxtype *vtxdist;	/* Distribution of vertices */

  idxtype *match;
  idxtype *cmap;

  /* Communication/Setup parameters */
  int nnbrs, nrecv, nsend;	/* The number of neighboring processors */
  idxtype *peind;		/* Array of size nnbrs storing the neighboring PEs */
  idxtype *sendptr, *sendind;	/* CSR format of the vertices that are sent */
  idxtype *recvptr, *recvind;	/* CSR format of the vertices that are received */
  idxtype *imap;		/* The inverse map of local to global indices */
  idxtype *pexadj, *peadjncy, 
          *peadjloc;		/* CSR format of the PEs each vertex is adjancent to */

  int nlocal;			/* Number of interior vertices */
  idxtype *lperm;		/* lperm[0:nlocal] points to interior vertices, the rest are interface */

  /* Communication parameters for projecting the partition. 
   * These are computed during CreateCoarseGraph and used during projection 
   * Note that during projection, the meaning of received and sent is reversed! */
  idxtype *rlens, *slens;	/* Arrays of size nnbrs of how many vertices you are sending and receiving */
  KeyValueType *rcand;


  /* Partition parameters */
  idxtype *where;		/* processor where vtx resides */
  idxtype *fusedinfo;		/* fused element vtx belongs to */
  idxtype *glblvtxid;
  idxtype *lpwgts, *gpwgts;
  realtype *lpvol, *gpvol;	/* The volume of the partitions */
  realtype *lpsurf, *gpsurf;	/* The surface area of the partitions */

  realtype lminratio, gminratio;	/* The value of the objective function */

  RInfoType *rinfo;

  /* Node refinement information */
  NRInfoType *nrinfo;

  realtype lmincut, mincut;
};

typedef struct graphdef MGridGraphType;


/*************************************************************************
* The following data type implements a timer
**************************************************************************/
typedef double timer;


/*************************************************************************
* The following structure stores information used by parallel kmetis
**************************************************************************/
struct controldef {
  int mype, npes;		/* Info about the parallel system */
  int dbglvl;			/* Controls the debuging output of the program */
  int nparts;			/* The number of partitions */

  int CType;                    /* The type of coarsening */
  int RType;                    /* The type of refinement */

  int minsize;                  /* The bounds on the number of elements per parti
tion */
  int maxsize;


  MPI_Comm gcomm;
  MPI_Comm comm;		/* MPI Communicator */
  MPI_Request sreq[MAX_PES], 
              rreq[MAX_PES];    /* MPI send and receive requests */
  MPI_Status statuses[MAX_PES];
  MPI_Status status;

  /* Various Timers */
  timer TotalTmr, InitPartTmr, MatchTmr, ContractTmr, CoarsenTmr, RefTmr,
        SetupTmr, ColorTmr, ProjectTmr, KWayInitTmr, KWayTmr, MoveTmr,
        RemapTmr, AuxTmr1, AuxTmr2, AuxTmr3, AuxTmr4, AuxTmr5, AuxTmr6;
};

typedef struct controldef MGridCtrlType;


/*************************************************************************
* The following data structure stores a sparse matrix in CSR format
* The diagonal entry is in the first position of each row.
**************************************************************************/
struct matrixdef {
  int nrows;		/* Number of rows in the matrix */
  idxtype *rowptr;
  idxtype *colind;
  float *values;
  idxtype *transfer;
};

typedef struct matrixdef MatrixType;
