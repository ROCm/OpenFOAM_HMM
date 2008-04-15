/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * refine.c
 *
 * This file contains the driving routines for multilevel k-way refinement
 *
 * George Irene
 */

#include "mgridgen.h"


/*************************************************************************
* This function is the entry point of refinement
**************************************************************************/
void RefineKWay(CtrlType *ctrl, GraphType *orggraph, GraphType *graph, int npasses)
{
  int i;

  ctrl->nparts = graph->nvtxs;

  graph->where = idxmalloc(graph->nvtxs, "graph->where");
  for (i=0; i<graph->nvtxs; i++)
     graph->where[i] = i;

  ComputeKWayPartitionParams(ctrl, graph);

  /* Get into the refinement loop */
  for (;;) {
     switch (ctrl->RType) {
       case REFINE_AR:
         Random_KWayARatioRefine(ctrl, graph, npasses); 
         break;
       case REFINE_WAR:
         Random_KWayWeightARatioRefine(ctrl, graph, npasses); 
         break;
       case REFINE_SCUT:
         Random_KWaySCutRefine(ctrl, graph, npasses); 
         break;
       case REFINE_MINMAXAVAR:
         Random_KWayMinMaxAverageARatioRefine(ctrl, graph, npasses); 
         break;
       case REFINE_MINMAXAR:
         Random_KWayMinMaxARatioRefine(ctrl, graph, npasses); 
         break;
       case REFINE_MULTIOBJECTIVE:
         Random_KWayMultiObjRefine(ctrl, graph, npasses); 
         break;
       case REFINE_MULTIOBJECTIVE2:
         Random_KWayMultiObjRefine2(ctrl, graph, npasses);
         break;
       default:
         errexit("Unknown RType of %d\n", ctrl->RType);
     }

     if (graph == orggraph) 
       break;
     else 
       graph = graph->finer;

     ProjectKWayPartition(graph);
     BreakComponents(ctrl, graph);

     Merge(ctrl, graph, npasses);

     ComputeKWayPartitionParams(ctrl, graph);

     IFSET(ctrl->dbglvl, DBG_REFINE,
           printf("Level done nparts=%d minratio=%e\n", ctrl->nparts,
                  graph->minratio));
  }

  BreakComponents(ctrl, graph);
  Merge(ctrl, graph, npasses);

  IMfree(&graph->pwgts, &graph->pvol, &graph->psurf, LTERM);
  ComputeKWayPartitionParams(ctrl, graph);
  Random_KWayMultiObjRefine(ctrl, graph, npasses); 

  Cycle(ctrl, graph, npasses);

  IMfree(&graph->pwgts, &graph->pvol, &graph->psurf, LTERM);
  IFSET(ctrl->dbglvl, DBG_REFINE, ComputeKWayPartitionParams(ctrl, graph));
  IFSET(ctrl->dbglvl, DBG_REFINE,
        printf("Last level done nparts=%d minratio=%e\n",
               ctrl->nparts, graph->minratio));

  IMfree(&graph->pwgts, &graph->pvol, &graph->psurf, LTERM);
  IFSET(ctrl->dbglvl, DBG_TRACK, ComputeKWayPartitionParams(ctrl, graph));
  IFSET(ctrl->dbglvl, DBG_TRACK, ComputeGridStatistics(ctrl, graph));
}


/*************************************************************************
* This function is the entry point of refinement
**************************************************************************/
void RefineKWayOnce(CtrlType *ctrl, GraphType *graph, int npasses)
{
  int i, nvtxs;
  int nparts, lastpart;
  idxKeyValueType *pairs;

  nvtxs = graph->nvtxs;

  pairs = (idxKeyValueType *) IMmalloc(nvtxs * sizeof(idxKeyValueType), "pairs");
  /* Find the number of elements: npart */
  for (i=0; i<nvtxs; i++) {
     pairs[i].key = graph->where[i];
     pairs[i].val = i;
  }
  idxkeysort(nvtxs, pairs);

  nparts = 1;
  lastpart = pairs[0].key;
  pairs[0].key = nparts - 1;
  for (i=1; i<nvtxs; i++) {
     if (pairs[i].key > lastpart) {
       nparts++;
       lastpart = pairs[i].key;
     }
     pairs[i].key = nparts - 1;
  }

  ctrl->nparts = nparts;

  for (i=0; i<nvtxs; i++)
     graph->where[pairs[i].val] = pairs[i].key;

  IMfree(&pairs, LTERM);

  /* Perform the refinement */
  ComputeKWayPartitionParams(ctrl, graph);

  switch (ctrl->RType) {
    case REFINE_AR:
      Random_KWayARatioRefine(ctrl, graph, npasses); 
      break;
    case REFINE_WAR:
      Random_KWayWeightARatioRefine(ctrl, graph, npasses); 
      break;
    case REFINE_SCUT:
      Random_KWaySCutRefine(ctrl, graph, npasses); 
      break;
    case REFINE_MINMAXAVAR:
      Random_KWayMinMaxAverageARatioRefine(ctrl, graph, npasses); 
      break;
    case REFINE_MINMAXAR:
      Random_KWayMinMaxARatioRefine(ctrl, graph, npasses); 
      break;
    case REFINE_MULTIOBJECTIVE:
      Random_KWayMultiObjRefine(ctrl, graph, npasses); 
      break;
    case REFINE_MULTIOBJECTIVE2:
      Random_KWayMultiObjRefine2(ctrl, graph, npasses); 
      break;
    default:
      errexit("Unknown RType of %d\n", ctrl->RType);
  }

  BreakComponents(ctrl, graph);
  Merge(ctrl, graph, npasses);

  IMfree(&graph->pwgts, &graph->pvol, &graph->psurf, LTERM);
  ComputeKWayPartitionParams(ctrl, graph);
  Random_KWayMultiObjRefine(ctrl, graph, npasses); 

  Cycle(ctrl, graph, npasses);

  IMfree(&graph->pwgts, &graph->pvol, &graph->psurf, LTERM);
  IFSET(ctrl->dbglvl, DBG_REFINE, ComputeKWayPartitionParams(ctrl, graph));
  IFSET(ctrl->dbglvl, DBG_REFINE,
        printf("Last level done nparts=%d minratio=%e\n", ctrl->nparts, graph->minratio));

  IMfree(&graph->pwgts, &graph->pvol, &graph->psurf, LTERM);
  IFSET(ctrl->dbglvl, DBG_TRACK, ComputeKWayPartitionParams(ctrl, graph));
  IFSET(ctrl->dbglvl, DBG_TRACK, ComputeGridStatistics(ctrl, graph));
}


/*************************************************************************
* This function computes the parameters required for the partitioning
**************************************************************************/
void ComputeKWayPartitionParams(CtrlType *ctrl, GraphType *graph)
{
  int i, j, me, nvtxs, nparts; 
  idxtype *xadj, *vwgt, *adjncy, *where, *pwgts;
  realtype *vvol, *vsurf, *adjwgt, *pvol, *psurf;


  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  vvol = graph->vvol;
  vsurf = graph->vsurf;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where = graph->where;

  nparts = ctrl->nparts;
  pwgts = graph->pwgts = idxsmalloc(nparts, 0, "pwgts");
  pvol = graph->pvol = realsmalloc(nparts, 0.0, "pvol");
  psurf = graph->psurf = realsmalloc(nparts, 0.0, "psurf");

  for (i=0; i<nvtxs; i++) {
     me = where[i];
     pwgts[me] += vwgt[i];
     pvol[me]  += vvol[i];
     psurf[me] += vsurf[i];

     for (j=xadj[i]; j<xadj[i+1]; j++) {
        if (where[adjncy[j]] != me)
          psurf[me] += adjwgt[j];
     }
  }

  graph->minratio = ComputeFunction(ctrl->RType, ctrl, graph);
}


/*************************************************************************
* This function projects a partition, and at the same time computes the
* parameters for refinement.
**************************************************************************/
void ProjectKWayPartition(GraphType *graph)
{
  int i, nvtxs;
  idxtype *where, *cwhere, *cmap;
  GraphType *cgraph;

  nvtxs = graph->nvtxs;
  cmap = graph->cmap;
  where = graph->where = idxmalloc(nvtxs, "where");

  cgraph = graph->coarser;
  cwhere = cgraph->where;

  /* Go through and project partition */
  for (i=0; i<nvtxs; i++) 
     where[i] = cwhere[cmap[i]];

  FreeGraph(cgraph);
  IMfree(&(graph->coarser), LTERM);
}

/*************************************************************************
* This function computes the statistics for the grid
**************************************************************************/
void ComputeGridStatistics(CtrlType *ctrl, GraphType *graph)
{
  int i, j, dim, nparts, nvtxs, from, to;
  idxtype *pwgts, *counts;
  idxtype *xadj, *adjncy, *where;
  realtype ed, min, max, sum, wsum, ratio, surf;
  realtype *pvol, *psurf;


  dim    = ctrl->dim;
  nparts = ctrl->nparts;
  where  = graph->where;
  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  adjncy = graph->adjncy;
  pwgts  = graph->pwgts;
  pvol   = graph->pvol;
  psurf  = graph->psurf;

  counts = idxsmalloc(ctrl->maxsize+1, 0, "counts");

  min = max = sum = ARATIO(dim, psurf[0], pvol[0]);
  wsum = 1.0*pwgts[0]*ARATIO(dim, psurf[0], pvol[0]);
  surf = psurf[0];
  counts[pwgts[0]]++;
  for (i=1; i<nparts; i++) {
     ratio = ARATIO(dim, psurf[i], pvol[i]);
     sum += ratio;
     wsum += 1.0*pwgts[i]*ratio;
     surf += psurf[i];
     if (min > ratio)
       min = ratio;
     if (max < ratio)
       max = ratio;
     counts[pwgts[i]]++;
  }

  for (ed=0.0, i=0; i<nvtxs; i++) {
     from = where[i];
     for (j=xadj[i]; j<xadj[i+1]; j++) {
        to = where[adjncy[j]];
        if (to != from)
          ed = ed + 1.0;
     }
  }
  ed = ed/2;

  printf("Npoints: %d, Coarsening Factor: %f\n", nparts, 1.0*graph->nvtxs/(1.0*nparts));
  printf("Aspect Ratios: Min : %e, Max : %e\n", min, max);
  printf("Aspect Ratios: Sum : %e, Wsum: %e\n", sum, wsum);
  printf("Aspect Ratios: Surf: %e, Avg : %e\n", surf, sum/(1.0*nparts));
  printf("Graph mincut : %e\n", ed);
  printf("Cell size: min=%d, max=%d\n", ctrl->minsize, ctrl->maxsize);
  printf("CellSizeDist: ");
  for (i=1; i<=ctrl->maxsize; i++)
     if (counts[i] != 0)
       printf("[%2d %4d] ", i, counts[i]);
  printf("\n");

  IMfree(&counts,LTERM);
}


/*************************************************************************
* This function finds all the connected components induced by the 
* partitioning vector and creates new partitions for each one of the 
* components
**************************************************************************/
void BreakComponents(CtrlType *ctrl, GraphType *graph)
{
  int i, j, k, me, nvtxs, nparts, first, last, nleft, ncmps;
  idxtype *xadj, *adjncy, *where;
  idxtype *touched, *perm, *todo, *cind, *cptr;

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  adjncy = graph->adjncy;
  where  = graph->where;

  nparts = ctrl->nparts;

  touched = idxsmalloc(nvtxs+1, 0, "touched");
  cptr    = idxmalloc(nvtxs+1, "cptr");
  cind    = idxmalloc(nvtxs+1, "cind");
  perm    = idxmalloc(nvtxs+1, "perm");
  todo    = idxmalloc(nvtxs+1, "todo");

  for (i=0; i<nvtxs; i++) 
     perm[i] = todo[i] = i;

  /* Find the connected componends induced by the partition */
  ncmps = -1;
  first = last = 0;
  nleft = nvtxs;
  while (nleft > 0) {
    if (first == last) { /* Find another starting vertex */
      cptr[++ncmps] = first;
      ASSERT(touched[todo[0]] == 0);
      i = todo[0];
      cind[last++] = i;
      touched[i] = 1;
      me = where[i];
    }

    i = cind[first++];
    k = perm[i];
    j = todo[k] = todo[--nleft];
    perm[j] = k;

    for (j=xadj[i]; j<xadj[i+1]; j++) {
       k = adjncy[j];
       if (where[k] == me && !touched[k]) {
         cind[last++] = k;
         touched[k] = 1;
       }
    }
  }
  cptr[++ncmps] = first;

  /* printf("I found %d components, for this %d-way partition\n", ncmps, nparts); */

  if (ncmps > nparts) {
    for (i=0; i<ncmps; i++) {
       for (j=cptr[i]; j<cptr[i+1]; j++)
          where[cind[j]] = i;
    }
    ctrl->nparts = ncmps;
  }

  IMfree(&touched, &cptr, &cind, &perm, &todo, LTERM);
}


/*************************************************************************
* This function computes the value of the RType function for the grid
**************************************************************************/
realtype ComputeFunction(int RType, CtrlType *ctrl, GraphType *graph)
{
  int i, dim, nparts;
  idxtype *pwgts;
  realtype new, ratio;
  realtype *pvol, *psurf;

  dim    = ctrl->dim;
  nparts = ctrl->nparts;
  pvol   = graph->pvol;
  psurf  = graph->psurf;
  pwgts  = graph->pwgts;

  switch (RType) {
    case REFINE_AR:
      ratio = ARATIO(dim, psurf[0], pvol[0]);
      for (i=1; i<nparts; i++)
         ratio += ARATIO(dim, psurf[i], pvol[i]);
      break;
    case REFINE_WAR:
      ratio = 1.0*pwgts[0]*ARATIO(dim, psurf[0], pvol[0]);
      for (i=1; i<nparts; i++)
         ratio += 1.0*pwgts[i]*ARATIO(dim, psurf[i], pvol[i]);
      break;
    case REFINE_SCUT:
      ratio = psurf[0];
      for (i=1; i<nparts; i++)
         ratio += psurf[i];
      break;
    case REFINE_MINMAXAVAR:
    case REFINE_MINMAXAR:
    case REFINE_MULTIOBJECTIVE:
    case REFINE_MULTIOBJECTIVE2:
      ratio = ARATIO(dim, psurf[0], pvol[0]);
      for (i=1; i<nparts; i++) {
         new = ARATIO(dim, psurf[i], pvol[i]);
        if (new > ratio)
          ratio = new;
      }
      break;
    default:
      errexit("Unknown RType of %d\n", ctrl->RType);
  }

  return(ratio);
}


/*************************************************************************
* This function computes the values of the RType functions for the grid
**************************************************************************/
void ComputeAllFunctions(CtrlType *ctrl, GraphType *graph)
{
  realtype ratio;

  ratio = ComputeFunction(1, ctrl, graph);
  printf("\t RTYPE=1 %e\n",ratio);
      
  ratio = ComputeFunction(2, ctrl, graph);
  printf("\t RTYPE=2 %e\n",ratio);
      
  ratio = ComputeFunction(3, ctrl, graph);
  printf("\t RTYPE=3 %e\n",ratio);
    
  ratio = ComputeFunction(4, ctrl, graph);
  printf("\t RTYPE=4 %e\n",ratio);

  ratio = ComputeFunction(5, ctrl, graph);
  printf("\t RTYPE=5 %e\n",ratio);
}
