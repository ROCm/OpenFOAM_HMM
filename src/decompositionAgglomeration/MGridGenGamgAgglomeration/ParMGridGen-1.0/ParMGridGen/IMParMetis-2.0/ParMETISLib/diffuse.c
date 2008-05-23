/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * diffuse.c
 *
 * This is the entry point of parallel difussive repartitioning routines
 *
 * George Irene
 */

#include "parmetis.h"


/***********************************************************************************
* This function is the entry point of the parallel multilevel local diffusion
* algorithm. It uses parallel undirected diffusion followed by adaptive k-way 
* refinement. This function utilizes local coarsening.
************************************************************************************/
void ParMETIS_RepartLDiffusion(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, 
       idxtype *vwgt, realtype *adjwgt, int *wgtflag, int *numflag, int *options,
       int *edgecut, idxtype *part, MPI_Comm *comm)
{
  int npes, mype;
  CtrlType ctrl;
  WorkSpaceType wspace;
  GraphType *graph;

  MPI_Comm_size(*comm, &npes);
  MPI_Comm_rank(*comm, &mype);

  if (npes == 1) { /* Take care the npes = 1 case */
    idxset(vtxdist[1], 0, part);
    *edgecut = 0;
    return;
  }

  if (*numflag == 1) 
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 1);

  SetUpCtrl(&ctrl, npes, options, *comm);
  ctrl.CoarsenTo = amin(vtxdist[npes]+1, 70*npes);

  graph = SetUpGraph(&ctrl, vtxdist, xadj, vwgt, adjncy, adjwgt, *wgtflag);
  graph->vsize = idxsmalloc(graph->nvtxs, 1, "Par_KMetis: vsize");

  PreAllocateMemory(&ctrl, graph, &wspace);

  IFSET(ctrl.dbglvl, DBG_TRACK, printf("%d ParMETIS_RepartLDiffusion about to call AdaptiveUndirected_Partition\n",mype));
  AdaptiveUndirected_Partition(&ctrl, graph, &wspace);

  IFSET(ctrl.dbglvl, DBG_TRACK, printf("%d ParMETIS_RepartLDiffusion about to call ReMapGraph\n",mype));
  ReMapGraph(&ctrl, graph, 0, &wspace);

  idxcopy(graph->nvtxs, graph->where, part);
  *edgecut = graph->mincut;

  IMfree(&graph->vsize, LTERM);
  FreeInitialGraphAndRemap(graph, *wgtflag);
  FreeWSpace(&wspace);
  FreeCtrl(&ctrl);

  if (*numflag == 1)
    ChangeNumbering(vtxdist, xadj, adjncy, part, npes, mype, 0);
}
