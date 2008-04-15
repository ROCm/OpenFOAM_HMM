/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * drivers.c
 *
 * This file contains the driving routines for the various parallel
 * multilevel partitioning and repartitioning algorithms
 *
 * George Irene
 */

#include "parmetis.h"


/*************************************************************************
* This function is the driver for the adaptive refinement mode of ParMETIS
**************************************************************************/
void AdaptiveUndirected_Partition(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace)
{
  int i;

  SetUp(ctrl, graph, wspace);

  IFSET(ctrl->dbglvl, DBG_PROGRESS, rprintf(ctrl, "[%6d %8d %5d %5d][%d][%d,%d]\n", 
        graph->gnvtxs, GlobalSESum(ctrl, graph->nedges), GlobalSEMin(ctrl, graph->nvtxs),
        GlobalSEMax(ctrl, graph->nvtxs), ctrl->CoarsenTo, graph->maxvwgt,
        GlobalSEMax(ctrl, graph->vwgt[idxamax(graph->nvtxs, graph->vwgt)])));  

  if (graph->gnvtxs < 1.3*ctrl->CoarsenTo || (graph->finer != NULL && graph->gnvtxs > graph->finer->gnvtxs*COARSEN_FRACTION)) {
    /* Set the initial partition */
    IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->InitPartTmr));
    graph->where = idxmalloc(graph->nvtxs+graph->nrecv, "Adaptive_Partition: graph->where");
    for (i=0; i<graph->nvtxs; i++)
      graph->where[i] = ctrl->mype;
    IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->InitPartTmr));

    ComputePartitionParams(ctrl, graph, wspace);
    KWayAdaptiveRefineClean(ctrl, graph, wspace, NGR_PASSES, UNBALANCE_FRACTION);

  }
  else { /* Coarsen it and the partition it */
    IFSET(ctrl->dbglvl, DBG_TRACK, printf("%d AdaptiveUndirected_Partition about to call LocalMatch_HEM\n",ctrl->mype));
    LocalMatch_HEM(ctrl, graph, wspace);

    AdaptiveUndirected_Partition(ctrl, graph->coarser, wspace);

    ProjectPartition(ctrl, graph, wspace);
    ComputePartitionParams(ctrl, graph, wspace);
    if (1.0*ctrl->nparts*graph->gpwgts[idxamax(ctrl->nparts, graph->gpwgts)]/(1.0*idxsum(ctrl->nparts, graph->gpwgts)) - UNBALANCE_FRACTION > 0.004)
      KWayAdaptiveRefineClean(ctrl, graph, wspace, NGR_PASSES, UNBALANCE_FRACTION);
    else
      KWayRefineClean(ctrl, graph, wspace, NGR_PASSES, UNBALANCE_FRACTION);
  }
}

