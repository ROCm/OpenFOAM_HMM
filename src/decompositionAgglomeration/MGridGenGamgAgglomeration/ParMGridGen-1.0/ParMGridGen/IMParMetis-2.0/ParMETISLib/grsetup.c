/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * grsetup.c
 *
 * This file contain various graph setting up routines
 *
 * George Irene
 */

#include "parmetis.h"



/*************************************************************************
* This function setsup the CtrlType structure
**************************************************************************/
GraphType *SetUpGraph(CtrlType *ctrl, idxtype *vtxdist, idxtype *xadj, 
              idxtype *vwgt, idxtype *adjncy, realtype *adjwgt, int wgtflag)
{
  int tvwgt;
  GraphType *graph;

  graph = CreateGraph();
  graph->level = 0;
  graph->gnvtxs = vtxdist[ctrl->npes];
  graph->nvtxs = vtxdist[ctrl->mype+1]-vtxdist[ctrl->mype];
  graph->nedges = xadj[graph->nvtxs];
  graph->xadj = xadj;
  graph->vwgt = vwgt;
  graph->adjncy = adjncy;
  graph->adjwgt = adjwgt;
  graph->vtxdist = vtxdist;

  if ((wgtflag&2) == 0)
    graph->vwgt = idxsmalloc(graph->nvtxs, 1, "Par_KMetis: vwgt");
  if ((wgtflag&1) == 0)
    graph->adjwgt = realsmalloc(graph->nedges, 1, "Par_KMetis: vwgt");

  tvwgt = GlobalSESum(ctrl, idxsum(graph->nvtxs, graph->vwgt));
  graph->maxvwgt = MAXVWGT_FACTOR*tvwgt/ctrl->CoarsenTo;

  return graph;
}


/*************************************************************************
* This function setsup the CtrlType structure
**************************************************************************/
void SetUpCtrl(CtrlType *ctrl, int nparts, int *options, MPI_Comm comm)
{
  MPI_Comm_dup(comm, &(ctrl->gcomm));
  MPI_Comm_rank(ctrl->gcomm, &ctrl->mype);
  MPI_Comm_size(ctrl->gcomm, &ctrl->npes);

  if (options[0] == 1) { /* Get user specified options */
    ctrl->dbglvl = options[OPTION_DBGLVL];
    ctrl->foldf = options[OPTION_FOLDF];
    ctrl->ipart = options[OPTION_IPART];
  }
  else { /* Assign default options */
    ctrl->dbglvl = 0;
    ctrl->foldf = 0;
    ctrl->ipart = IPART_RB;
  }

  ctrl->nparts = nparts;    /* Set the # of partitions is de-coupled from the # of domains */
  ctrl->comm = ctrl->gcomm;
  ctrl->xyztype = XYZ_SPFILL;

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
* This function computes movement statistics for adaptive refinement
* schemes
**************************************************************************/
void ComputeMoveStatistics(CtrlType *ctrl, GraphType *graph, int *nmoved, int *maxin, int *maxout)
{
  int i, j, nvtxs;
  idxtype *where;
  idxtype *lpvtxs, *gpvtxs;

  nvtxs = graph->nvtxs;
  where = graph->where;

  lpvtxs = idxsmalloc(ctrl->nparts, 0, "ComputeMoveStatistics: lpvtxs");
  gpvtxs = idxsmalloc(ctrl->nparts, 0, "ComputeMoveStatistics: gpvtxs");

  for (j=i=0; i<nvtxs; i++) {
    lpvtxs[where[i]]++;
    if (where[i] != ctrl->mype)
      j++;
  }

  /* PrintVector(ctrl, ctrl->npes, 0, lpvtxs, "Lpvtxs: "); */

  MPI_Allreduce((void *)lpvtxs, (void *)gpvtxs, ctrl->nparts, IDX_DATATYPE, MPI_SUM, ctrl->comm);

  *nmoved = GlobalSESum(ctrl, j);
  *maxout = GlobalSEMax(ctrl, j);
  *maxin = GlobalSEMax(ctrl, gpvtxs[ctrl->mype]-(nvtxs-j));

  IMfree(&lpvtxs, &gpvtxs, LTERM);
}
