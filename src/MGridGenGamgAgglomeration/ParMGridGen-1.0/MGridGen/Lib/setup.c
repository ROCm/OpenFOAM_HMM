/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * setup.c
 *
 * This file contains setup functions for mgridgen
 *
 * George Irene
 */

#include "mgridgen.h"


/*************************************************************************
* This function sets up the graph from the user input
**************************************************************************/
void SetUpGraph(GraphType *graph, int nvtxs, idxtype *xadj, realtype *vvol,
                realtype *vsurf, idxtype *adjncy, realtype *adjwgt)
{
  int i, j;

  graph->nvtxs     = nvtxs;
  graph->xadj      = idxmalloc(nvtxs+1, "xadj");
  graph->vwgt      = idxsmalloc(nvtxs, 1, "vwgt");
  graph->vvol      = realmalloc(nvtxs, "vvol");
  graph->vsurf     = realmalloc(nvtxs, "vsurf");
  graph->adjncy    = idxmalloc(xadj[nvtxs], "adjncy");
  graph->adjwgt    = realmalloc(xadj[nvtxs], "adjwgt");
  graph->adjwgtsum = realsmalloc(nvtxs, 0.0, "adjwgtsum");

  graph->pwgts = NULL;
  graph->pvol = NULL;
  graph->psurf = NULL;

  idxcopy(nvtxs+1, xadj, graph->xadj);
  realcopy(nvtxs, vvol, graph->vvol);
  realcopy(nvtxs, vsurf, graph->vsurf);
  idxcopy(xadj[nvtxs], adjncy, graph->adjncy);
  realcopy(xadj[nvtxs], adjwgt, graph->adjwgt);

  for (i=0; i<nvtxs; i++)
     for (j=xadj[i]; j<xadj[i+1]; j++)
        graph->adjwgtsum[i] += adjwgt[j];
}


/*************************************************************************
* This function frees the memory that was allocated for the graph
**************************************************************************/
void FreeGraph(GraphType *graph)
{
  IMfree(&graph->xadj, &graph->vwgt, &graph->vvol, &graph->vsurf,
              &graph->adjncy, &graph->adjwgt, &graph->adjwgtsum, &graph->cmap,
              &graph->where, &graph->pwgts, &graph->pvol, &graph->psurf, LTERM);

}


/*************************************************************************
* This function creates a CoarseGraphType data structure and initializes
* the various fields
**************************************************************************/
GraphType *CreateGraph(void)
{
  GraphType *graph;

  graph = (GraphType *)IMmalloc(sizeof(GraphType), "CreateGraph: graph");

  graph->nvtxs = -1;
  graph->nmoves = -1;
  graph->xadj = NULL;
  graph->vwgt = NULL;
  graph->vvol = NULL;
  graph->vsurf = NULL;
  graph->adjncy = NULL;
  graph->adjwgt = NULL;
  graph->adjwgtsum = NULL;
  graph->cmap = NULL;
  graph->where = NULL;
  graph->pwgts = NULL;
  graph->pvol = NULL;
  graph->psurf = NULL;

  return graph;
}
