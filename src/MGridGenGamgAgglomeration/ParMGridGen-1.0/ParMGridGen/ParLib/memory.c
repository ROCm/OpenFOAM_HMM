/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * memory.c
 *
 * This file contains routines that deal with memory allocation
 *
 * George Irene
 */

#include "parmgridgen.h"


/*************************************************************************
* This function allocate various pools of memory
**************************************************************************/
void PreAllocateMGridMemory(MGridCtrlType *ctrl, MGridGraphType *graph, MGridWorkSpaceType *wspace)
{
  wspace->nlarge  = 2*graph->nedges;

/* pairs  : 2*nedges. KeyValueType: 2*idxtype.   So 4*nedges(idxtype) */
/* indices: 2*nedges.                            So 2*nedges(idxtype) */
/* degrees: 1*nedges. EdgeType: idxtype+realtype.So 1*nedges(idxtype)+1*nedges(realtype) */

  wspace->maxcore = graph->nedges*(2*sizeof(KeyValueType) + 2*sizeof(idxtype) + sizeof(EdgeType));

  wspace->core = IMmalloc(wspace->maxcore, "PreAllocateMemory: wspace->core");

  wspace->pairs = (KeyValueType *)wspace->core;
  wspace->indices = (idxtype *)(wspace->pairs + wspace->nlarge);
  wspace->degrees = (EdgeType *)(wspace->indices + wspace->nlarge);


  wspace->pv1 = idxmalloc(ctrl->npes+1, "PreAllocateMemory: wspace->pv1");
  wspace->pv2 = idxmalloc(ctrl->npes+1, "PreAllocateMemory: wspace->pv2");
  wspace->pv4 = idxmalloc(ctrl->npes+1, "PreAllocateMemory: wspace->pv4");

  wspace->pepairs1 = (KeyValueType *)IMmalloc(sizeof(KeyValueType)*(ctrl->npes+1), "PreAllocateMemory: wspace->pepairs1");
  wspace->pepairs2 = (KeyValueType *)IMmalloc(sizeof(KeyValueType)*(ctrl->npes+1), "PreAllocateMemory: wspace->pepairs2");

}


/*************************************************************************
* This function de-allocate various pools of memory
**************************************************************************/
void FreeMGridWSpace(MGridWorkSpaceType *wspace)
{

  IMfree(&wspace->core, 
         &wspace->pv1, 
         &wspace->pv2, 
         &wspace->pv4, 
         &wspace->pepairs1, 
         &wspace->pepairs2, 
         LTERM);
}


/*************************************************************************
* This function de-allocates memory allocated for the control structures
**************************************************************************/
void FreeMGridCtrl(MGridCtrlType *ctrl)
{
  MPI_Comm_free(&(ctrl->gcomm));
}


/*************************************************************************
* This function creates a MGridGraphType data structure and initializes
* the various fields
**************************************************************************/
MGridGraphType *CreateMGridGraph(void)
{
  MGridGraphType *graph;

  graph = (MGridGraphType *)IMmalloc(sizeof(MGridGraphType), "CreateMGridGraph: graph");

  InitMGridGraph(graph);

  return graph;
}


/*************************************************************************
* This function creates a CoarseGraphType data structure and initializes
* the various fields
**************************************************************************/
void InitMGridGraph(MGridGraphType *graph) 
{
  graph->gnvtxs = graph->nvtxs = graph->nedges = -1;
  graph->nnbrs = graph->nrecv = graph->nsend = graph->nlocal = -1;
  graph->xadj = graph->vwgt = graph->vsize = graph->adjncy = NULL;
  graph->vvol = graph->vsurf = NULL;
  graph->lpvol = graph->gpvol = graph->lpsurf = graph->gpsurf = NULL;
  graph->adjwgt = graph->adjwgtsum = NULL;
  graph->lminratio = graph->gminratio = -1.0;
  graph->vtxdist = NULL;
  graph->match = graph->cmap = NULL;

  graph->peind = NULL;
  graph->sendptr = graph->sendind = graph->recvptr = graph->recvind = NULL;
  graph->imap = NULL;
  graph->pexadj = graph->peadjncy = graph->peadjloc = NULL;
  graph->lperm = NULL;

  graph->slens = graph->rlens = NULL;
  graph->rcand = NULL;

  graph->glblvtxid = graph->fusedinfo = graph->where = graph->lpwgts = graph->gpwgts = NULL;
  graph->rinfo = NULL;

  graph->nrinfo = NULL;
}

/*************************************************************************
* This function deallocates any memory stored in a graph
**************************************************************************/
void FreeMGridGraph(MGridGraphType *graph) 
{

  IMfree(&graph->xadj, 
         &graph->vwgt,
         &graph->vvol,
         &graph->vsurf,
         &graph->adjwgtsum,
         &graph->lpvol,
         &graph->gpvol,
         &graph->lpsurf,
         &graph->gpsurf,
         &graph->vsize,
         &graph->adjncy,
         &graph->adjwgt,
         &graph->vtxdist, 
         &graph->match, 
         &graph->cmap, 
         &graph->lperm, 
         &graph->glblvtxid,
         &graph->fusedinfo, 
         &graph->where, 
         &graph->rinfo, 
         &graph->nrinfo, 
         &graph->lpwgts, 
         &graph->gpwgts, 
         &graph->peind, 
         &graph->sendptr, 
         &graph->sendind, 
         &graph->recvptr, 
         &graph->recvind, 
         &graph->imap,
         &graph->rlens,
         &graph->slens,
         &graph->rcand,
         &graph->pexadj,
         &graph->peadjncy,
         &graph->peadjloc,
         LTERM);

  IMfree(&graph, LTERM);
}


/*************************************************************************
* This function deallocates any memory stored in a graph
**************************************************************************/
void FreeMGridGraphContent(MGridGraphType *graph) 
{

  IMfree(&graph->xadj, 
         &graph->vwgt,
         &graph->vvol,
         &graph->vsurf,
         &graph->adjwgtsum,
         &graph->lpvol,
         &graph->gpvol,
         &graph->lpsurf,
         &graph->gpsurf,
         &graph->vsize,
         &graph->adjncy,
         &graph->adjwgt,
         &graph->vtxdist, 
         &graph->match, 
         &graph->cmap, 
         &graph->lperm, 
         &graph->glblvtxid,
         &graph->fusedinfo, 
         &graph->where, 
         &graph->rinfo, 
         &graph->nrinfo, 
         &graph->lpwgts, 
         &graph->gpwgts, 
         &graph->peind, 
         &graph->sendptr, 
         &graph->sendind, 
         &graph->recvptr, 
         &graph->recvind, 
         &graph->imap,
         &graph->rlens,
         &graph->slens,
         &graph->rcand,
         &graph->pexadj,
         &graph->peadjncy,
         &graph->peadjloc,
         LTERM);
}
