/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * coarsen.c
 *
 * This file contains the driving routines for the coarsening process 
 *
 * George Irene
 */

#include "mgridgen.h"


/*************************************************************************
* This function takes a graph and creates a sequence of coarser graphs
**************************************************************************/
GraphType *Coarsen(CtrlType *ctrl, GraphType *graph)
{
  int j=0;
  GraphType *cgraph;

  cgraph = graph;

  do {
    IFSET(ctrl->dbglvl, DBG_COARSEN, printf("%6d %7d\n", cgraph->nvtxs,
          cgraph->xadj[cgraph->nvtxs]));

    switch (ctrl->CType) {
      case MATCH_RM:
        Match_RM(ctrl, cgraph);
        break;
      case MATCH_HEM:
        Match_HEM(ctrl, cgraph);
        break;
      case MATCH_HEM_SLOW:
        Match_HEM_Slow(ctrl, cgraph);
        break;
      case MATCH_HEM_TRUE:
        Match_HEM_True(ctrl, cgraph);
        break;
      default:
        errexit("Unknown CType: %d\n", ctrl->CType);
    }

    j++;
    cgraph = cgraph->coarser;

  } while (cgraph->nvtxs < cgraph->finer->nvtxs);

  IFSET(ctrl->dbglvl, DBG_COARSEN, printf("Coarsening Info : %d %d %d\n",j,
        cgraph->nvtxs, cgraph->finer->nvtxs));
  return cgraph;
}


/*************************************************************************
* This function takes a graph and creates a sequence of coarser graphs
**************************************************************************/
GraphType *Coarsen_Restricted(CtrlType *ctrl, GraphType *graph)
{
  int i, nvtxs;
  idxtype *cmap, *where;
  GraphType *cgraph;

  cgraph = graph;

  do {
    IFSET(ctrl->dbglvl, DBG_COARSEN, printf("%6d %7d\n", cgraph->nvtxs,
    cgraph->xadj[cgraph->nvtxs]));
    Match_HEM_Slow_Restricted(ctrl, cgraph);

    /* Propagate the where vector downwards */
    nvtxs = cgraph->nvtxs;
    cmap  = cgraph->cmap;
    where = cgraph->where;

    cgraph = cgraph->coarser;
    cgraph->where = idxmalloc(cgraph->nvtxs, "cgraph->where");

    for (i=0; i<nvtxs; i++)
       cgraph->where[cmap[i]] = where[i];

    IMfree(&(cgraph->finer->where), LTERM);

  } while (cgraph->nvtxs < cgraph->finer->nvtxs);


  if (cgraph->nvtxs != ctrl->nparts) 
    printf("It appears that some domains are non-contigous [%d %d]\n",
           cgraph->nvtxs, ctrl->nparts);

  IMfree(&(cgraph->where), LTERM);


  /* Perform any additional coarsening */
  do {
    IFSET(ctrl->dbglvl, DBG_COARSEN, printf("%6d %7d\n", cgraph->nvtxs,
          cgraph->xadj[cgraph->nvtxs]));
    Match_HEM_Slow(ctrl, cgraph);

    cgraph = cgraph->coarser;

  } while (cgraph->nvtxs < cgraph->finer->nvtxs);

  return cgraph;
}


/*************************************************************************
* This function creates the coarser graph
**************************************************************************/
void CreateCoarseGraph(GraphType *graph, int cnvtxs, idxtype *match,
                       idxtype *perm)
{
  int i, j, k, m, nvtxs, nedges, cnedges, v, u;
  idxtype *xadj, *vwgt, *adjncy;
  idxtype *cxadj, *cvwgt, *cadjncy;
  realtype *vvol, *vsurf, *adjwgt, *adjwgtsum;
  realtype *cvvol, *cvsurf, *cadjwgt, *cadjwgtsum;
  idxtype *cmap, *htable;
  GraphType *cgraph;


  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  vvol = graph->vvol;
  vsurf = graph->vsurf;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  adjwgtsum = graph->adjwgtsum;
  cmap = graph->cmap;

  /* Initialize the coarser graph */
  cgraph = SetUpCoarseGraph(graph, cnvtxs);
  cxadj      = cgraph->xadj;
  cvwgt      = cgraph->vwgt;
  cvvol      = cgraph->vvol;
  cvsurf     = cgraph->vsurf;
  cadjncy    = cgraph->adjncy;
  cadjwgt    = cgraph->adjwgt;
  cadjwgtsum = cgraph->adjwgtsum;

  htable = idxsmalloc(cnvtxs, -1, "htable");

  cxadj[0] = cnvtxs = cnedges = 0;
  for (i=0; i<nvtxs; i++) {
     v = perm[i];
     if (cmap[v] != cnvtxs) 
       continue;

     u = match[v];
     cvwgt[cnvtxs]      = vwgt[v];
     cvvol[cnvtxs]      = vvol[v];
     cvsurf[cnvtxs]     = vsurf[v];
     cadjwgtsum[cnvtxs] = adjwgtsum[v];

     for (nedges=0, j=xadj[v]; j<xadj[v+1]; j++) {
        k = cmap[adjncy[j]];
        if ((m = htable[k]) == -1) {
          cadjncy[nedges] = k;
          cadjwgt[nedges] = adjwgt[j];
          htable[k] = nedges++;
        }
        else
          cadjwgt[m] += adjwgt[j];
     }

     if (v != u) { 
       cvwgt[cnvtxs]      += vwgt[u];
       cvvol[cnvtxs]      += vvol[u];
       cvsurf[cnvtxs]     += vsurf[u];
       cadjwgtsum[cnvtxs] += adjwgtsum[u];

       for (j=xadj[u]; j<xadj[u+1]; j++) {
          k = cmap[adjncy[j]];
          if ((m = htable[k]) == -1) {
            cadjncy[nedges] = k;
            cadjwgt[nedges] = adjwgt[j];
            htable[k] = nedges++;
          }
          else
            cadjwgt[m] += adjwgt[j];
       }

       /* Remove the contracted adjacency weight */
       if ((j = htable[cnvtxs]) != -1) {
         ASSERT(cadjncy[j] == cnvtxs);
         cadjwgtsum[cnvtxs] -= cadjwgt[j];
         cadjncy[j] = cadjncy[--nedges];
         cadjwgt[j] = cadjwgt[nedges];
         htable[cnvtxs] = -1;
       }
     }

     for (j=0; j<nedges; j++)
        htable[cadjncy[j]] = -1;  /* Zero out the htable */

     cnedges += nedges;
     cxadj[++cnvtxs] = cnedges;
     cadjncy += nedges;
     cadjwgt += nedges;
  }

  free(htable);
}


/*************************************************************************
* Setup the various arrays for the coarse graph
**************************************************************************/
GraphType *SetUpCoarseGraph(GraphType *graph, int cnvtxs)
{
  int nedges;
  GraphType *cgraph;

  cgraph = CreateGraph();
  cgraph->nvtxs = cnvtxs;
  cgraph->finer = graph;
  graph->coarser = cgraph;

  nedges = graph->xadj[graph->nvtxs];

  /* Allocate memory for the coarser graph */
  cgraph->xadj       = idxmalloc(cnvtxs+1, "xadj");
  cgraph->vwgt       = idxmalloc(cnvtxs, "vwgt");
  cgraph->vvol       = realmalloc(cnvtxs, "vvol");
  cgraph->vsurf      = realmalloc(cnvtxs, "vsurf");
  cgraph->adjwgtsum  = realmalloc(cnvtxs, "adjwgtsum");
  cgraph->adjncy     = idxmalloc(nedges, "adjncy");
  cgraph->adjwgt     = realmalloc(nedges, "adjwgt");


  return cgraph;
}
