/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * move.c
 *
 * This file contains functions that move the graph given a partition
 *
 * George Irene
 */

#include "parmgridgen.h"

/*************************************************************************
* This function moves the graph, and returns a new graph.
* This routine can be called with or without performing refinement.
* In the latter case it allocates and computes lpwgts itself.
**************************************************************************/
MGridGraphType *MoveMGridGraph(MGridCtrlType *ctrl, MGridGraphType *graph, MGridWorkSpaceType *wspace)
{
  int i, ii, j, jj, nvtxs, nparts, maxrealcore;
  idxtype *xadj, *vwgt, *adjncy, *mvtxdist, *p_d;
  idxtype *where, *newlabel, *lpwgts, *gpwgts, *fusedinfo, *glblvtxid;
  idxtype *sidxgraph, *ridxgraph;
  realtype *adjwgt, *vvol, *vsurf;
  realtype *srealgraph, *rrealgraph, *realcore;
  KeyValueType *sinfo, *rinfo;
  MGridGraphType *mgraph;

  nparts = ctrl->nparts;
  ASSERT(ctrl, nparts == ctrl->npes);

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  vvol = graph->vvol;
  vsurf = graph->vsurf;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where = graph->where;
  fusedinfo = graph->fusedinfo;
  glblvtxid = graph->glblvtxid;

  mvtxdist = idxmalloc(nparts+1, "MoveMGridGraph: mvtxdist");

  /* Let's do a prefix scan to determine the labeling of the nodes given */
  lpwgts = wspace->pv1;
  gpwgts = wspace->pv2;
  sinfo = wspace->pepairs1;
  rinfo = wspace->pepairs2;
  for (i=0; i<nparts; i++)
    sinfo[i].key = sinfo[i].val = 0;

  /* Here we care about the count and not total weight (diff since graph may be weighted) */
  for (i=0; i<nvtxs; i++) {
    ASSERTP(ctrl, where[i] >= 0 && where[i] < ctrl->npes, (ctrl, "%d %d %d\n", i, where[i], ctrl->npes) );
    sinfo[where[i]].key++;
    sinfo[where[i]].val += xadj[i+1]-xadj[i];
  }
  for (i=0; i<nparts; i++)
    lpwgts[i] = sinfo[i].key;

  MPI_Scan((void *)lpwgts, (void *)gpwgts, nparts, IDX_DATATYPE, MPI_SUM, ctrl->comm);
  MPI_Allreduce((void *)lpwgts, (void *)mvtxdist, nparts, IDX_DATATYPE, MPI_SUM, ctrl->comm);

  MAKECSR(i, nparts, mvtxdist);

  /* gpwgts[i] will store the label of the first vertex for each domain in each processor */
  for (i=0; i<nparts; i++)
    gpwgts[i] = mvtxdist[i] + gpwgts[i] - lpwgts[i];  /* We were interested in an exclusive Scan */

  newlabel = idxmalloc(nvtxs+graph->nrecv, "MoveMGridGraph: newlabel");

  for (i=0; i<nvtxs; i++) 
    newlabel[i] = gpwgts[where[i]]++;

  /* OK, now send the newlabel info to processors storing adjacent interface nodes */
  MGridCommInterfaceData(ctrl, graph, newlabel, wspace->indices, newlabel+nvtxs);

  /* Now lets tell everybody what and from where he will get it. Assume nparts == npes */
  MPI_Alltoall((void *)sinfo, 2, IDX_DATATYPE, (void *)rinfo, 2, IDX_DATATYPE, ctrl->comm);

  /* Use lpwgts and gpwgts as pointers to where data will be received and sent */
  lpwgts[0] = 0;  /* Send part */
  gpwgts[0] = 0;  /* Received part */
  for (i=0; i<nparts; i++) {
    lpwgts[i+1] = lpwgts[i] + 3*sinfo[i].key + sinfo[i].val;
    gpwgts[i+1] = gpwgts[i] + 3*rinfo[i].key + rinfo[i].val;
  }


  if (lpwgts[nparts]+gpwgts[nparts] > wspace->maxcore) {
    /* Adjust core memory, incase the graph was originally very memory unbalanced */
    free(wspace->core);
    wspace->maxcore = lpwgts[nparts]+4*gpwgts[nparts]; /* In spirit of the 8*nedges */
    wspace->core = idxmalloc(wspace->maxcore, "MoveMGridGraph: wspace->core");
  }

  maxrealcore = (1 + sizeof(realtype)/sizeof(idxtype))*graph->nedges;
  if (lpwgts[nparts]+gpwgts[nparts] > maxrealcore) {
    maxrealcore = lpwgts[nparts]+4*gpwgts[nparts]; /* In spirit of the 8*nedges */
    realcore = realmalloc(maxrealcore, "MoveMGridGraph: realcore");
  }
  else
    realcore = realmalloc(maxrealcore, "PreAllocateMemory: realcore");

  sidxgraph = wspace->core;
  ridxgraph = wspace->core + lpwgts[nparts];

  srealgraph = realcore;
  rrealgraph = realcore + lpwgts[nparts];

  /* Issue the receives first */
  for (i=0; i<nparts; i++) {
    if (rinfo[i].key > 0) 
      MPI_Irecv((void *)(ridxgraph+gpwgts[i]), gpwgts[i+1]-gpwgts[i], IDX_DATATYPE, i, 1, ctrl->comm, ctrl->rreq+i);
    else
      ASSERT(ctrl, gpwgts[i+1]-gpwgts[i] == 0);
  }

  /* Assemble the graph to be sent and send it */
  for (i=0; i<nvtxs; i++) {
    ii = lpwgts[where[i]];
    sidxgraph[ii] = xadj[i+1]-xadj[i];
    srealgraph[ii++] = vvol[i];
    sidxgraph[ii] = vwgt[i]; 
    srealgraph[ii++] = vsurf[i];
    sidxgraph[ii] = fusedinfo[i]; 
/*    srealgraph[ii++] = -1;  */
    p_d = (idxtype *) &srealgraph[ii++];    /* srealgraph[ii++] = glblvtxid[i] is NOT safe ! */
    *p_d = glblvtxid[i];
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      sidxgraph[ii] = newlabel[adjncy[j]];
      srealgraph[ii++] = adjwgt[j];
    }
    lpwgts[where[i]] = ii;
  }

  for (i=nparts; i>0; i--)
    lpwgts[i] = lpwgts[i-1];
  lpwgts[0] = 0;

  for (i=0; i<nparts; i++) {
    if (sinfo[i].key > 0)
      MPI_Isend((void *)(sidxgraph+lpwgts[i]), lpwgts[i+1]-lpwgts[i], IDX_DATATYPE, i, 1, ctrl->comm, ctrl->sreq+i);
    else
      ASSERT(ctrl, lpwgts[i+1]-lpwgts[i] == 0);
  }

/*
#ifdef DMALLOC
  ASSERT(ctrl, dmalloc_verify(NULL) == DMALLOC_VERIFY_NOERROR);
#endif
*/

  /* Wait for the send/recv to finish */
  for (i=0; i<nparts; i++) {
    if (sinfo[i].key > 0) 
      MPI_Wait(ctrl->sreq+i, &ctrl->status);
  }
  for (i=0; i<nparts; i++) {
    if (rinfo[i].key > 0) 
      MPI_Wait(ctrl->rreq+i, &ctrl->status); 
  }


  /* Issue the receives first */
  for (i=0; i<nparts; i++) {
    if (rinfo[i].key > 0) 
      MPI_Irecv((void *)(rrealgraph+gpwgts[i]), gpwgts[i+1]-gpwgts[i], REAL_DATATYPE, i, 1, ctrl->comm, ctrl->rreq+i);
    else
      ASSERT(ctrl, gpwgts[i+1]-gpwgts[i] == 0);
  }

  for (i=0; i<nparts; i++) {
    if (sinfo[i].key > 0)
      MPI_Isend((void *)(srealgraph+lpwgts[i]), lpwgts[i+1]-lpwgts[i], REAL_DATATYPE, i, 1, ctrl->comm, ctrl->sreq+i);
    else
      ASSERT(ctrl, lpwgts[i+1]-lpwgts[i] == 0);
  }

/*
#ifdef DMALLOC
  ASSERT(ctrl, dmalloc_verify(NULL) == DMALLOC_VERIFY_NOERROR);
#endif
*/

  /* Wait for the send/recv to finish */
  for (i=0; i<nparts; i++) {
    if (sinfo[i].key > 0) 
      MPI_Wait(ctrl->sreq+i, &ctrl->status);
  }
  for (i=0; i<nparts; i++) {
    if (rinfo[i].key > 0) 
      MPI_Wait(ctrl->rreq+i, &ctrl->status); 
  }

  /* OK, now go and put the graph into GraphType Format */
  mgraph = CreateMGridGraph();
  mgraph->gnvtxs = graph->gnvtxs;
  mgraph->nvtxs = mgraph->nedges = 0;
  for (i=0; i<nparts; i++) {
    mgraph->nvtxs += rinfo[i].key;
    mgraph->nedges += rinfo[i].val;
  }
  nvtxs = mgraph->nvtxs;
  xadj = mgraph->xadj = idxmalloc(nvtxs+1, "MoveGraph: mgraph->xadj");
  vwgt = mgraph->vwgt = idxmalloc(nvtxs, "MoveGraph: mgraph->vwgt");
  vvol = mgraph->vvol = realmalloc(nvtxs, "MoveGraph: mgraph->vvol");
  vsurf = mgraph->vsurf = realmalloc(nvtxs, "MoveGraph: mgraph->vsurf");
  mgraph->adjwgtsum = realmalloc(nvtxs, "MoveGraph: mgraph->adjwgtsum");
  adjncy = mgraph->adjncy = idxmalloc(mgraph->nedges, "MoveGraph: mgraph->adjncy");
  adjwgt = mgraph->adjwgt = realmalloc(mgraph->nedges, "MoveGraph: mgraph->adjwgt");
  fusedinfo = mgraph->fusedinfo = idxmalloc(nvtxs, "MoveGraph: mgraph->fusedinfo");
  glblvtxid = mgraph->glblvtxid = idxmalloc(nvtxs, "MoveGraph: mgraph->glblvtxid");
  mgraph->vtxdist = mvtxdist;

  for (jj=ii=i=0; i<nvtxs; i++) {
    xadj[i] = ridxgraph[ii];
    vvol[i] = rrealgraph[ii++];
    vwgt[i] = ridxgraph[ii]; 
    vsurf[i] = rrealgraph[ii++];
    fusedinfo[i] = ridxgraph[ii];
    p_d = (idxtype *) &rrealgraph[ii++];
    glblvtxid[i] = *p_d;
    for (j=0; j<xadj[i]; j++) {
      adjncy[jj] = ridxgraph[ii];
      adjwgt[jj++] = rrealgraph[ii++];
    }
  }
  MAKECSR(i, nvtxs, xadj);

  for (i=0; i<mgraph->nvtxs; i++) {
    mgraph->adjwgtsum[i] = 0.0;
    for (j=xadj[i]; j<xadj[i+1]; j++)
      mgraph->adjwgtsum[i] += adjwgt[j];
  }

  ASSERTP(ctrl, jj == mgraph->nedges, (ctrl, "%d %d\n", jj, mgraph->nedges));
  ASSERTP(ctrl, ii == gpwgts[nparts], (ctrl, "%d %d %d %d %d\n", ii, gpwgts[nparts], jj, mgraph->nedges, nvtxs));

  IMfree(&newlabel, &realcore, LTERM);
  /* CheckMGraph(ctrl, mgraph); */

  return mgraph;
}

