/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * grsetup.c
 *
 * This file contain various graph setting up routines
 *
 * George Irene
 */

#include "parmgridgen.h"


/*************************************************************************
* This function setsup the GraphType structure
**************************************************************************/
MGridGraphType *SetUpMGridGraph(MGridCtrlType *ctrl, idxtype *vtxdist, idxtype *xadj, realtype *vvol,
                      realtype *vsurf, idxtype *adjncy, realtype *adjwgt)
{
  int i, j;
  MGridGraphType *graph;

  graph = CreateMGridGraph();
  graph->gnvtxs = vtxdist[ctrl->npes];
  graph->nvtxs = vtxdist[ctrl->mype+1]-vtxdist[ctrl->mype];
  graph->nedges = xadj[graph->nvtxs];

  graph->vtxdist = vtxdist;
  graph->xadj = xadj;
  graph->vvol = vvol;
  graph->vsurf = vsurf;
  graph->adjncy = adjncy;
  graph->adjwgt = adjwgt;

  graph->vwgt = idxsmalloc(graph->nvtxs, 1, "SetUpMGridGraph: vwgt");

  graph->adjwgtsum = realsmalloc(graph->nvtxs, 0.0, "adjwgtsum");
  for (i=0; i<graph->nvtxs; i++)
    for (j=xadj[i]; j<xadj[i+1]; j++)
      graph->adjwgtsum[i] += adjwgt[j];

/*  tvwgt = MGridGlobalSESum(ctrl, MGrididxsum(graph->nvtxs, graph->vwgt)); */

  return graph;
}


/*************************************************************************
* This function setsup the CtrlType structure
**************************************************************************/
void SetUpMGridCtrl(MGridCtrlType *ctrl, int minsize, int maxsize, int *options, MPI_Comm comm)
{
  MPI_Comm_dup(comm, &(ctrl->gcomm));
  MPI_Comm_rank(ctrl->gcomm, &ctrl->mype);
  MPI_Comm_size(ctrl->gcomm, &ctrl->npes);

  ctrl->minsize = minsize;
  ctrl->maxsize = maxsize;
  ctrl->CType = options[OPTION_CTYPE];
  ctrl->RType = options[OPTION_RTYPE];
  ctrl->dbglvl = options[OPTION_DBGLVL];

  ctrl->comm = ctrl->gcomm;

  srand(ctrl->mype);
  srand48(ctrl->mype);
}


/*************************************************************************
* This function changes the numbering from 1 to 0 or 0 to 1
**************************************************************************/
void ChangeNumbering(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *part, int npes, int mype, int from)
{
  int i, nvtxs, nedges;

  if (from == 1) {  /* Change it from 1 to 0 */
    for (i=0; i<npes+1; i++)
      vtxdist[i]--;

    nvtxs = vtxdist[mype+1]-vtxdist[mype];
    for (i=0; i<nvtxs+1; i++) 
      xadj[i]--;

    nedges = xadj[nvtxs];
    for (i=0; i<nedges; i++) 
      adjncy[i]--;
  }
  else {  /* Change it from 0 to 1 */
    nvtxs = vtxdist[mype+1]-vtxdist[mype];
    nedges = xadj[nvtxs];

    for (i=0; i<npes+1; i++) 
      vtxdist[i]++;

    for (i=0; i<nvtxs+1; i++) 
      xadj[i]++; 

    for (i=0; i<nedges; i++) 
      adjncy[i]++; 

    for (i=0; i<nvtxs; i++)
      part[i]++;

  }
}


/*************************************************************************
* This function removes outside edges and makes graph purely local
**************************************************************************/
void RemoveEdges(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy,
                 realtype *adjwgt, idxtype *newxadj, idxtype *newadjncy,
                 realtype *newadjwgt, MPI_Comm comm)
{
  int  i, j, k, l, mype;
  int firstvtx, lastvtx, nvtxs;

  MPI_Comm_rank(comm, &mype);

  firstvtx = vtxdist[mype];
  lastvtx = vtxdist[mype+1];
  nvtxs = lastvtx - firstvtx;

  for (newxadj[0]=0, l=0, i=0; i<nvtxs; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (k>=firstvtx && k<lastvtx) {
        newadjncy[l] = k-firstvtx;
        newadjwgt[l++] = adjwgt[j];
      }
    }
    newxadj[i+1]=l;
  }

}

/*************************************************************************
* This function corrects surface edges (vsurf)
**************************************************************************/
void CorrectSurfaceEdges(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy,
                         realtype *adjwgt, realtype *vsurf, MPI_Comm comm)
{
  int  i, j, k, mype;
  int firstvtx, lastvtx, nvtxs;

  MPI_Comm_rank(comm, &mype);

  firstvtx = vtxdist[mype];
  lastvtx = vtxdist[mype+1];
  nvtxs = lastvtx - firstvtx;

  for (i=0; i<nvtxs; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (k<firstvtx || k>=lastvtx)
        vsurf[i] += adjwgt[j];
    }
  }

}

/*************************************************************************
* This function removes outside edges and corrects surface edges (vsurf)
**************************************************************************/
void RemoveCorrectSurfaceEdges(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy,
               realtype *adjwgt, idxtype *newxadj, idxtype *newadjncy,
               realtype *newadjwgt, realtype *newvsurf, MPI_Comm comm)
{
  int  i, j, k, l, mype;
  int firstvtx, lastvtx, nvtxs;

  MPI_Comm_rank(comm, &mype);

  firstvtx = vtxdist[mype];
  lastvtx = vtxdist[mype+1];
  nvtxs = lastvtx - firstvtx;

  for (newxadj[0]=0, l=0, i=0; i<nvtxs; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];
      if (k>=firstvtx && k<lastvtx) {
        newadjncy[l] = k-firstvtx;
        newadjwgt[l++] = adjwgt[j];
      }
      else
        newvsurf[i] += adjwgt[j];
    }
    newxadj[i+1]=l;
  }

}
