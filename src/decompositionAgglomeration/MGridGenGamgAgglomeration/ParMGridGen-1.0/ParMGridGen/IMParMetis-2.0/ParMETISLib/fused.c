/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * fused.c
 *
 * This is the creation of the fused element graph
 *
 * Started 04/18/00
 * Irene
 */

#include "parmetis.h"

/***********************************************************************************
* This function creates the fused-element-graph and returns the partition
************************************************************************************/
void ParMETIS_FusedElementGraph(idxtype *vtxdist, idxtype *xadj, realtype *vvol,
              realtype *vsurf, idxtype *adjncy, idxtype *vwgt, realtype *adjwgt,
              int *wgtflag, int *numflag, int *nparts, int *options,
              idxtype *part, MPI_Comm *comm)
{
  int npes, mype, nvtxs;
  CtrlType ctrl;
  WorkSpaceType wspace;
  GraphType *graph;

  MPI_Comm_size(*comm, &npes);
  MPI_Comm_rank(*comm, &mype);
 
  nvtxs = vtxdist[mype+1]-vtxdist[mype];

  /* IFSET(options[OPTION_DBGLVL], DBG_TRACK, printf("%d ParMETIS_FEG npes=%d\n",mype, npes)); */

  SetUpCtrl(&ctrl, *nparts, options, *comm);
  ctrl.CoarsenTo = amin(vtxdist[npes]+1, 25*amax(npes, *nparts));

  graph = SetUpGraph(&ctrl, vtxdist, xadj, vwgt, adjncy, adjwgt, *wgtflag);

  graph->where = part;

  PreAllocateMemory(&ctrl, graph, &wspace);

  IFSET(ctrl.dbglvl, DBG_TIME, InitTimers(&ctrl));
  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, starttimer(ctrl.TotalTmr));

  CreateFusedElementGraph(&ctrl, graph, &wspace, numflag);

  idxcopy(nvtxs, graph->where, part);

  IFSET(ctrl.dbglvl, DBG_TIME, MPI_Barrier(ctrl.gcomm));
  IFSET(ctrl.dbglvl, DBG_TIME, stoptimer(ctrl.TotalTmr));

  if (((*wgtflag)&2) == 0)
    IMfree(&graph->vwgt, LTERM);
  IMfree(&graph->lperm, &graph->peind, &graph->pexadj, &graph->peadjncy,
         &graph->peadjloc, &graph->recvptr, &graph->recvind, &graph->sendptr,
         &graph->imap, &graph->sendind, &graph, LTERM);
  FreeWSpace(&wspace);
  FreeCtrl(&ctrl);
}

/**********************************************************************************/
void CreateFusedElementGraph(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace, int *numflag)
{
  int i, j, k, l;
  int *nparts, ipart, mypart, newpart;
  int wgtflag, edgecut, npes, mype, nvtxs, nedges, counter;
  int gnfvtxs, nfvtxs, firstfvtx;
  int foptions[10];
  idxtype *fptr, *find;
  idxtype *part, *fpart, *spart, *rpart;
  idxtype *vtxdist, *xadj, *adjncy;
  idxtype *fvtxdist, *fxadj, *fadjncy;
  idxtype *map;
  realtype *fadjwgt;
  MPI_Comm *comm;

  npes = ctrl->npes;
  mype = ctrl->mype;
  nparts = &(ctrl->nparts);
  comm = &(ctrl->comm);

  vtxdist = graph->vtxdist;
  xadj = graph->xadj;
  adjncy = graph->adjncy;

  nvtxs = vtxdist[mype+1]-vtxdist[mype];
  nedges = xadj[nvtxs];

  SetUp(ctrl, graph, wspace);

  /* Communicate number of parts found */
  fvtxdist = idxmalloc(npes+1, "FusedElementGraph: fvtxdist");
  MPI_Allgather((void *)nparts, 1, MPI_INT, (void *)fvtxdist, 1, MPI_INT, *comm);

  MAKECSR(i, npes, fvtxdist); 
  firstfvtx = fvtxdist[mype];
  nfvtxs = fvtxdist[mype+1]-fvtxdist[mype];
  gnfvtxs = fvtxdist[npes];

  ASSERT(ctrl, nfvtxs == *nparts);

  part = idxmalloc(nvtxs+graph->nrecv, "FusedElementGraph: part");
  idxcopy(nvtxs, graph->where, part);
  spart = wspace->indices;
  rpart = part + nvtxs;

  CommInterfaceData(ctrl, graph, part, spart, rpart);

  /* Create a part-to-vertex mapping */
/*  map = idxsmalloc(gnfvtxs, -1, "FusedElementGraph: map"); */ /* TOO GENEROUS !! */
/*  fptr = idxsmalloc(*nparts+1, 0, "FusedElementGraph: fptr"); */
/*  find = idxmalloc(nvtxs, "FusedElementGraph: find"); */
  if (gnfvtxs + nvtxs + (*nparts+1) <= wspace->maxcore) {
    map = wspace->core;
    idxset(gnfvtxs, -1, map);
    fptr = map + gnfvtxs;
  }
  else {
    map = idxsmalloc(gnfvtxs, -1, "FusedElementGraph: map"); 
    fptr = wspace->core;
  }
  idxset((*nparts+1), 0, fptr);
  find = fptr + (*nparts+1);

  for (i=0; i<nvtxs; i++)
     fptr[part[i]-firstfvtx]++;
  MAKECSR(i, *nparts, fptr);
  for (i=0; i<nvtxs; i++) {
     ipart = part[i] - firstfvtx;
     find[fptr[ipart]] = i;
     fptr[ipart]++;
  }

  for (ipart=*nparts; ipart>0; ipart--)
     fptr[ipart] = fptr[ipart-1];
  fptr[0] = 0;

  /* Create the fused graph for the local edges */
  fxadj = idxsmalloc(nfvtxs+1, 0, "FusedElementGraph: fxadj");
  fadjncy = idxmalloc(nedges, "FusedElementGraph: fadjncy");
  fadjwgt = realsmalloc(nedges, 0, "FusedElementGraph: fadjwgt");

  fxadj[0] = 0;
  for (ipart=0; ipart<*nparts; ipart++) {
     counter = 0;
     mypart = ipart + firstfvtx;

     for (l=fptr[ipart]; l<fptr[ipart+1]; l++) {
        i = find[l];
        for (j=xadj[i]; j<xadj[i+1]; j++) {
           k=adjncy[j];
           newpart=part[k];
           if (newpart != mypart && map[newpart] == -1) {  /* edge is not created yet */
             map[newpart] = fxadj[ipart]+counter; 
             fadjncy[fxadj[ipart]+counter] = newpart;
             fadjwgt[map[newpart]] = 1;         /* alternatively = adjwgt[k] */
             counter++;
           }
           else if (newpart != mypart && map[newpart] != -1)  /* edge is already there */
             fadjwgt[map[newpart]]++;
        }
     }
     fxadj[ipart+1] = fxadj[ipart] + counter;
     for (i=fxadj[ipart]; i<fxadj[ipart+1]; i++)
        map[fadjncy[i]] = -1;
  }

  /* Now change the weights of the interface edges */
  ChangeWeights(nfvtxs, fvtxdist, fxadj, fadjncy, fadjwgt, *comm);
  
  /* Repartition the graph using fused elements */
  foptions[0] = 1;
  foptions[3] = 0;
  
/* fpart = idxmalloc(nfvtxs, "TestParMetis: fpart"); */
  fpart = map;   /* it is OK since nfvtxs < gnfvtxs */
  
  wgtflag = 1;
  ParMETIS_RepartLDiffusion(fvtxdist, fxadj, fadjncy, NULL, fadjwgt,
                 &wgtflag, numflag, foptions, &edgecut, fpart, comm);
  
  /* Project the partitioning back to the original graph */
  for (ipart=0; ipart<nfvtxs; ipart++)  {
     ASSERTP(ctrl, fpart[ipart] >= 0 && fpart[ipart] < npes, (ctrl, "%d %d %d\n", ipart , fpart[ipart], npes) );
     for (i=fptr[ipart]; i<fptr[ipart+1]; i++)
        graph->where[find[i]]=fpart[ipart];       
  }

  if (gnfvtxs + nvtxs + (*nparts+1) > wspace->maxcore)
    IMfree(&map, LTERM);
  IMfree(&fvtxdist, &fxadj, &fadjncy, &fadjwgt, &part, LTERM);
}
