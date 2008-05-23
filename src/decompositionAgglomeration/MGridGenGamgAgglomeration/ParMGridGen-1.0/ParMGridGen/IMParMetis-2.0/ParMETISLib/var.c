/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * var.c
 *
 * This file contains routines related to armetis
 *
 * Irene
 */

#include "parmetis.h"

/*************************************************************************
* This function changes the weights of interface edges
**************************************************************************/
void ChangeWeights(int nvtxs, idxtype *vtxdist, idxtype *xadj, idxtype *adjncy,
                 realtype *adjwgt, MPI_Comm comm)
{
  int  i, j, k, mype;
  int firstvtx, lastvtx;

  MPI_Comm_rank(comm, &mype);

  firstvtx = vtxdist[mype];
  lastvtx = vtxdist[mype+1];

  for (i=0; i<nvtxs; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (k<firstvtx || k>=lastvtx)
        adjwgt[j] *= 10;       /* PLAY HERE !! */
    }
  }

}
