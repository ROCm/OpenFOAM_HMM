/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * refine.c
 *
 * This file contains code that performs the k-way refinement
 *
 * George Irene
 */

#include "parmetis.h"

#define DEBUG_PROJECT_
#define DEBUG_COMPUTEPPARAM_


#define ProperSide(c, from, other) \
              (((c) == 0 && (from)-(other) < 0) || ((c) == 1 && (from)-(other) > 0))

/*************************************************************************
* This function projects a partition.
**************************************************************************/
void ProjectPartition(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace)
{
  int i, nvtxs, nnbrs, firstvtx, cfirstvtx;
  idxtype *match, *cmap, *where, *cwhere;
  idxtype *peind, *slens, *rlens;
  KeyValueType *rcand, *scand;
  GraphType *cgraph;


  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->ProjectTmr));

  cgraph = graph->coarser;
  cwhere = cgraph->where;
  cfirstvtx = cgraph->vtxdist[ctrl->mype];

  nvtxs = graph->nvtxs;
  match = graph->match;
  cmap = graph->cmap;
  where = graph->where = idxmalloc(nvtxs+graph->nrecv, "ProjectPartition: graph->where");
  firstvtx = graph->vtxdist[ctrl->mype];


  if (graph->match_type == MATCH_GLOBAL) {  /* Only if global matching is on */
    /*------------------------------------------------------------
    / Start the transmission of the remote where information 
    /------------------------------------------------------------*/
    scand = wspace->pairs;
    nnbrs = graph->nnbrs;
    peind = graph->peind;
    slens = graph->slens;
    rlens = graph->rlens;
    rcand = graph->rcand;

    /* Issue the receives first */
    for (i=0; i<nnbrs; i++) {
      if (slens[i+1]-slens[i] > 0) /* Issue a receive only if you are getting something */
        MPI_Irecv((void *)(scand+slens[i]), 2*(slens[i+1]-slens[i]), IDX_DATATYPE, peind[i], 1, ctrl->comm, ctrl->rreq+i);
    }

#ifdef DEBUG_PROJECT
    PrintPairs(ctrl, rlens[nnbrs], rcand, "rcand"); 
#endif

    /* Put the where[rcand[].key] into the val field */
    for (i=0; i<rlens[nnbrs]; i++) {
      ASSERT(ctrl, rcand[i].val >= 0 && rcand[i].val < cgraph->nvtxs);
      rcand[i].val = cwhere[rcand[i].val];
    }

#ifdef DEBUG_PROJECT
    PrintPairs(ctrl, rlens[nnbrs], rcand, "rcand");
    PrintVector(ctrl, nvtxs, firstvtx, cmap, "cmap");
#endif

    /* Issue the sends next */
    for (i=0; i<nnbrs; i++) {
      if (rlens[i+1]-rlens[i] > 0) /* Issue a send only if you are sending something */
        MPI_Isend((void *)(rcand+rlens[i]), 2*(rlens[i+1]-rlens[i]), IDX_DATATYPE, peind[i], 1, ctrl->comm, ctrl->sreq+i);
    }
  }

  /*------------------------------------------------------------
  / Project local vertices first
  /------------------------------------------------------------*/
  for (i=0; i<nvtxs; i++) {
    if (match[i] >= KEEP_BIT) {
      ASSERT(ctrl, cmap[i]-cfirstvtx>=0 && cmap[i]-cfirstvtx<cgraph->nvtxs);
      where[i] = cwhere[cmap[i]-cfirstvtx];
    }
  }

  if (graph->match_type == MATCH_GLOBAL) {  /* Only if global matching is on */
    /*------------------------------------------------------------
    / Wait for the nonblocking operations to finish
    /------------------------------------------------------------*/
    for (i=0; i<nnbrs; i++) {
      if (rlens[i+1]-rlens[i] > 0)  
        MPI_Wait(ctrl->sreq+i, &ctrl->status);
    }
    for (i=0; i<nnbrs; i++) {
      if (slens[i+1]-slens[i] > 0)  
        MPI_Wait(ctrl->rreq+i, &ctrl->status);
    }

#ifdef DEBUG_PROJECT
    PrintPairs(ctrl, slens[nnbrs], scand, "scand"); 
#endif

    /*------------------------------------------------------------
    / Project received vertices now
    /------------------------------------------------------------*/
    for (i=0; i<slens[nnbrs]; i++) {
      ASSERTP(ctrl, scand[i].key-firstvtx>=0 && scand[i].key-firstvtx<graph->nvtxs, (ctrl, "%d %d %d\n", scand[i].key, firstvtx, graph->nvtxs));
      where[scand[i].key-firstvtx] = scand[i].val;
    }
  }


  FreeGraph(graph->coarser);
  graph->coarser = NULL;

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->ProjectTmr));
}



/*************************************************************************
* This function computes the initial id/ed 
**************************************************************************/
void ComputePartitionParams(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace)
{
  int i, j, k, nvtxs, firstvtx;
  idxtype *xadj, *ladjncy, *vtxdist, *lpwgts, *gpwgts;
  realtype *adjwgt;
  idxtype *where, *swhere, *rwhere;
  RInfoType *rinfo, *myrinfo;
  EdgeType *edegrees;
  int me, other;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->KWayInitTmr));

  nvtxs = graph->nvtxs;

  vtxdist = graph->vtxdist;
  xadj = graph->xadj;
  ladjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  where = graph->where;
  rinfo = graph->rinfo = (RInfoType *)IMmalloc(sizeof(RInfoType)*nvtxs, "ComputePartitionParams: rinfo");
  lpwgts = graph->lpwgts = idxsmalloc(ctrl->nparts, 0, "ComputePartitionParams: lpwgts");
  gpwgts = graph->gpwgts = idxmalloc(ctrl->nparts, "ComputePartitionParams: gpwgts");

  firstvtx = vtxdist[ctrl->mype];

  /*------------------------------------------------------------
  / Send/Receive the where information of interface vertices
  /------------------------------------------------------------*/
  swhere = wspace->indices;
  rwhere = where + nvtxs;

  CommInterfaceData(ctrl, graph, where, swhere, rwhere); 

#ifdef DEBUG_COMPUTEPPARAM
  PrintVector(ctrl, nvtxs, firstvtx, where, "where");
#endif

  ASSERT(ctrl, wspace->nlarge >= xadj[nvtxs]);

  /*------------------------------------------------------------
  / Compute now the id/ed degrees
  /------------------------------------------------------------*/
  graph->lmincut = 0;
  for (i=0; i<nvtxs; i++) {
    me = where[i];
    myrinfo = rinfo+i;

    lpwgts[me] += graph->vwgt[i];

    myrinfo->degrees = wspace->degrees + xadj[i];
    myrinfo->ndegrees = 0;
    myrinfo->id = myrinfo->ed = 0.0;

    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (me == where[ladjncy[j]])
        myrinfo->id += adjwgt[j];
      else
        myrinfo->ed += adjwgt[j];
    }


    if (myrinfo->ed > 0.0) {  /* Time to do some serious work */
      graph->lmincut += myrinfo->ed;
      edegrees = myrinfo->degrees;

      for (j=xadj[i]; j<xadj[i+1]; j++) {
        other = where[ladjncy[j]];
        if (me != other) {
          for (k=0; k<myrinfo->ndegrees; k++) {
            if (edegrees[k].edge == other) {
              edegrees[k].ewgt += adjwgt[j];
              break;
            }
          }
          if (k == myrinfo->ndegrees) {
            edegrees[k].edge = other;
            edegrees[k].ewgt = adjwgt[j];
            myrinfo->ndegrees++;
          }
          ASSERT(ctrl, myrinfo->ndegrees <= xadj[i+1]-xadj[i]);
        }
      }
    }
  }

#ifdef DEBUG_COMPUTEPPARAM
  PrintVector(ctrl, ctrl->nparts, 0, lpwgts, "lpwgts");
#endif

  /* Finally, sum-up the partition weights */
  MPI_Allreduce((void *)lpwgts, (void *)gpwgts, ctrl->nparts, IDX_DATATYPE, MPI_SUM, ctrl->comm);

  graph->mincut = GlobalSESumReal(ctrl, graph->lmincut)/2;

#ifdef DEBUG_COMPUTEPPARAM
  PrintVector(ctrl, ctrl->nparts, 0, gpwgts, "gpwgts");
#endif

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->KWayInitTmr));
}


/*************************************************************************
* This function performs k-way refinement
**************************************************************************/
void KWayRefineClean(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace, int npasses, float ubfraction)
{
  int i, ii, j, k, pass, nvtxs, firstvtx, otherlastvtx, c, nmoves, 
      nlupd, nsupd, nnbrs, nchanged;
  int npes = ctrl->npes, mype = ctrl->mype, nparts = ctrl->nparts;
  idxtype *xadj, *ladjncy, *vtxdist;
  realtype *adjwgt;
  idxtype *where, *lpwgts, *gpwgts;
  idxtype *peind, *recvptr, *sendptr;
  idxtype *update, *supdate, *rupdate, *pe_updates, *htable;
  idxtype *changed, *pperm, *perm;
  KeyValueType *swchanges, *rwchanges;
  int *nupds_pe;
  RInfoType *rinfo, *myrinfo;
  EdgeType *degrees;
  int from, me, other, vwgt, badminpwgt, badmaxpwgt;
  realtype oldcut;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->KWayTmr));

  nvtxs = graph->nvtxs;

  vtxdist = graph->vtxdist;
  xadj = graph->xadj;
  ladjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  firstvtx = vtxdist[mype];

  where = graph->where;
  rinfo = graph->rinfo;
  lpwgts = graph->lpwgts;
  gpwgts = graph->gpwgts;

  nnbrs = graph->nnbrs;
  peind = graph->peind;
  recvptr = graph->recvptr;
  sendptr = graph->sendptr;

  changed = idxmalloc(nvtxs, "KWayRefine: changed");
  rwchanges = wspace->pairs;
  swchanges = rwchanges + recvptr[nnbrs];

  update = idxmalloc(nvtxs, "KWayRefine: update");
  supdate = wspace->indices;
  rupdate = supdate + recvptr[nnbrs];
  nupds_pe = imalloc(npes, "KWayRefine: nupds_pe");

  htable = idxsmalloc(nvtxs+graph->nrecv, 0, "KWayRefine: lhtable");

  badminpwgt = (1.0/ubfraction)*idxsum(nparts, gpwgts)/nparts;
  badmaxpwgt = ubfraction*idxsum(nparts, gpwgts)/nparts;

  perm = idxmalloc(nvtxs, "KWayRefine: perm");
  for (i=0; i<nvtxs; i++)
    perm[i] = i;
  pperm = idxmalloc(nparts, "KWayRefine: pperm");
  for (i=0; i<nparts; i++)
    pperm[i] = i;

  IFSET(ctrl->dbglvl, DBG_REFINEINFO, rprintf(ctrl, "K-way refinement %2d [%5d, %5d], [%5d, %5d] Cut: %f, \tBalance: %6.3f\n",
        graph->level, gpwgts[idxamin(nparts, gpwgts)], gpwgts[idxamax(nparts, gpwgts)], badminpwgt, badmaxpwgt, graph->mincut,
        1.0*nparts*gpwgts[idxamax(nparts, gpwgts)]/(1.0*idxsum(nparts, gpwgts))));

  for (pass=0; pass<npasses; pass++) {
    if (mype == 0)
      RandomPermute(nparts, pperm, 0);
    MPI_Bcast((void *)pperm, nparts, IDX_DATATYPE, 0, ctrl->comm);

    oldcut = graph->mincut;

    /* PrintVector(ctrl, nvtxs, firstvtx, where, "where"); */

    FastRandomPermute(nvtxs, perm, 0);
    for (c=0; c<2; c++) {
      nlupd = nsupd = nmoves = nchanged = 0;
      for (ii=0; ii<nvtxs; ii++) {
        i = perm[ii];
        if (rinfo[i].ed >= rinfo[i].id) { /* Total ED is too high */
          degrees = rinfo[i].degrees;
          from = where[i];
          vwgt = graph->vwgt[i];

          if (gpwgts[from]-vwgt < badminpwgt)
            continue;   /* This cannot be moved! */

          for (k=0; k<rinfo[i].ndegrees; k++) {
            other = degrees[k].edge;
            if (ProperSide(c, pperm[from], pperm[other]) && gpwgts[other]+vwgt <= badmaxpwgt)
              break;
          }

          if (k < rinfo[i].ndegrees) { /* You actually found one */
            for (j=k+1; j<rinfo[i].ndegrees; j++) {
              other = degrees[j].edge;
              if (ProperSide(c, pperm[from], pperm[other]) && 
                  ((degrees[j].ewgt > degrees[k].ewgt && gpwgts[other]+vwgt <= badmaxpwgt) ||
                  (degrees[j].ewgt == degrees[k].ewgt && gpwgts[other] < gpwgts[degrees[k].edge])))
                k = j;
            }

            other = degrees[k].edge;

            if (degrees[k].ewgt > rinfo[i].id || 
                /* (degrees[k].ewgt-rinfo[i].id > -vwgt && gpwgts[other]+vwgt < 0.8*badminpwgt) || */
                (degrees[k].ewgt == rinfo[i].id && 
                 (other == mype || gpwgts[from] > badmaxpwgt || gpwgts[other] < badminpwgt || (from != mype && (gpwgts[from] - gpwgts[other] >= vwgt))))) {
              /* Update where, weight, and ID/ED information of the vertex you moved */
              where[i] = other;
              INC_DEC(lpwgts[other], lpwgts[from], vwgt);
              INC_DEC(gpwgts[other], gpwgts[from], vwgt);

              if (htable[i] == 0) {  /* make sure you do the update */
                htable[i] = 1;
                update[nlupd++] = i;
              }

              /* Put the vertices adjacent to i into the update array */
              for (j=xadj[i]; j<xadj[i+1]; j++) {
                k = ladjncy[j];
                if (htable[k] == 0) {
                  htable[k] = 1;
                  if (k<nvtxs)
                    update[nlupd++] = k;
                  else
                    supdate[nsupd++] = k;
                }
              }
              nmoves++;
              if (graph->pexadj[i+1]-graph->pexadj[i] > 0)
                changed[nchanged++] = i;
            }
          }
        }
      }

      /* myprintf(ctrl, "nmoves: %d, nlupd: %d, nsupd: %d\n", nmoves, nlupd, nsupd); */

      /* Tell everybody interested what the new where[] info is for the interface vertices */
      CommChangedInterfaceData(ctrl, graph, nchanged, changed, where, swchanges, rwchanges, wspace->pv4); 


      IFSET(ctrl->dbglvl, DBG_RMOVEINFO, rprintf(ctrl, "\t[%d %d], [%5d, %5d],  [%d %d %d]\n", 
                pass, c, badminpwgt, badmaxpwgt, 
                GlobalSESum(ctrl, nmoves), GlobalSESum(ctrl, nsupd), GlobalSESum(ctrl, nlupd)));


      /*-------------------------------------------------------------
      / Time to communicate with processors to send the vertices
      / whose degrees need to be update.
      /-------------------------------------------------------------*/
      /* Issue the receives first */
      for (i=0; i<nnbrs; i++) {
        MPI_Irecv((void *)(rupdate+sendptr[i]), sendptr[i+1]-sendptr[i], IDX_DATATYPE,
                  peind[i], 1, ctrl->comm, ctrl->rreq+i);
      }

      /* Issue the sends next. This needs some preporcessing */
      for (i=0; i<nsupd; i++) {
        htable[supdate[i]] = 0;
        supdate[i] = graph->imap[supdate[i]];
      }
      iidxsort(nsupd, supdate);

      for (j=i=0; i<nnbrs; i++) {
        otherlastvtx = vtxdist[peind[i]+1];
        for (k=j; k<nsupd && supdate[k] < otherlastvtx; k++); 
        MPI_Isend((void *)(supdate+j), k-j, IDX_DATATYPE, peind[i], 1, ctrl->comm, ctrl->sreq+i);
        j = k;
      }

      /* OK, now get into the loop waiting for the send/recv operations to finish */
      MPI_Waitall(nnbrs, ctrl->rreq, ctrl->statuses);
      for (i=0; i<nnbrs; i++) 
        MPI_Get_count(ctrl->statuses+i, IDX_DATATYPE, nupds_pe+i);
      MPI_Waitall(nnbrs, ctrl->sreq, ctrl->statuses);


      /*-------------------------------------------------------------
      / Place the recieved to-be updated vertices into update[] 
      /-------------------------------------------------------------*/
      for (i=0; i<nnbrs; i++) {
        pe_updates = rupdate+sendptr[i];
        for (j=0; j<nupds_pe[i]; j++) {
          k = pe_updates[j];
          if (htable[k-firstvtx] == 0) {
            htable[k-firstvtx] = 1;
            update[nlupd++] = k-firstvtx;
          }
        }
      }


      /*-------------------------------------------------------------
      / Update the rinfo of the vertices in the update[] array
      /-------------------------------------------------------------*/
      for (ii=0; ii<nlupd; ii++) {
        i = update[ii];
        ASSERT(ctrl, htable[i] == 1);

        htable[i] = 0;

        me = where[i];
        myrinfo = rinfo+i;

        graph->lmincut -= myrinfo->ed;

        myrinfo->ndegrees = 0;
        myrinfo->id = myrinfo->ed = 0.0;
        degrees = myrinfo->degrees;

        for (j=xadj[i]; j<xadj[i+1]; j++) {
          other = where[ladjncy[j]];
          if (me != other) {
            myrinfo->ed += adjwgt[j];

            for (k=0; k<myrinfo->ndegrees; k++) {
              if (degrees[k].edge == other) {
                degrees[k].ewgt += adjwgt[j];
                break;
              }
            }
            if (k == myrinfo->ndegrees) {
              degrees[k].edge = other;
              degrees[k].ewgt = adjwgt[j];
              myrinfo->ndegrees++;
            }
            ASSERT(ctrl, myrinfo->ndegrees <= xadj[i+1]-xadj[i]);
          }
          else {
            myrinfo->id += adjwgt[j];
          }
        }
        
        graph->lmincut += myrinfo->ed;
      }
    }

    /* Finally, sum-up the partition weights */
    MPI_Allreduce((void *)lpwgts, (void *)gpwgts, nparts, IDX_DATATYPE, MPI_SUM, ctrl->comm);

    graph->mincut = GlobalSESumReal(ctrl, graph->lmincut)/2;


    IFSET(ctrl->dbglvl, DBG_REFINEINFO, rprintf(ctrl, "\t[%5d, %5d], [%5d, %5d] Cut: %f\n", 
          gpwgts[idxamin(nparts, gpwgts)], gpwgts[idxamax(nparts, gpwgts)], badminpwgt, badmaxpwgt, graph->mincut));

    if (graph->mincut == oldcut && ubfraction > 1.00)
      break;
  }

  IMfree(&update, &nupds_pe, &htable, &changed, &pperm, &perm, LTERM);

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->KWayTmr));
}


/*************************************************************************
* This function performs k-way refinement
**************************************************************************/
void KWayAdaptiveRefineClean(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace, int npasses, float ubfraction)
{
  int i, ii, j, k, pass, nvtxs, firstvtx, otherlastvtx, c, nmoves, 
      nlupd, nsupd, nnbrs, nchanged;
  int npes = ctrl->npes, mype = ctrl->mype, nparts = ctrl->nparts;
  idxtype *xadj, *ladjncy, *vtxdist;
  realtype *adjwgt;
  idxtype *where, *lpwgts, *gpwgts;
  idxtype *peind, *recvptr, *sendptr;
  idxtype *update, *supdate, *rupdate, *pe_updates, *htable;
  idxtype *changed, *pperm, *perm, *tocheck;
  KeyValueType *swchanges, *rwchanges;
  int *nupds_pe;
  RInfoType *rinfo, *myrinfo;
  EdgeType *degrees;
  int from, me, other, to, vwgt, vsize, badminpwgt, badmaxpwgt, ntodo, ndone;
  realtype oldcut;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->KWayTmr));

  nvtxs = graph->nvtxs;

  vtxdist = graph->vtxdist;
  xadj = graph->xadj;
  ladjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  firstvtx = vtxdist[mype];

  where = graph->where;
  rinfo = graph->rinfo;
  lpwgts = graph->lpwgts;
  gpwgts = graph->gpwgts;

  nnbrs = graph->nnbrs;
  peind = graph->peind;
  recvptr = graph->recvptr;
  sendptr = graph->sendptr;

  changed = idxmalloc(nvtxs, "KWayRefine: changed");
  rwchanges = wspace->pairs;
  swchanges = rwchanges + recvptr[nnbrs];

  update = idxmalloc(nvtxs, "KWayRefine: update");
  supdate = wspace->indices;
  rupdate = supdate + recvptr[nnbrs];
  nupds_pe = imalloc(npes, "KWayRefine: nupds_pe");

  htable = idxsmalloc(nvtxs+graph->nrecv, 0, "KWayRefine: lhtable");

  tocheck = idxmalloc(nvtxs, "KWayRefine: tocheck");

  badminpwgt = (1.0/ubfraction)*idxsum(nparts, gpwgts)/nparts;
  badmaxpwgt = ubfraction*idxsum(nparts, gpwgts)/nparts;

  perm = idxmalloc(nvtxs, "KWayRefine: perm");
  for (i=0; i<nvtxs; i++)
    perm[i] = i;
  pperm = idxmalloc(nparts, "KWayRefine: pperm");
  for (i=0; i<nparts; i++)
    pperm[i] = i;

  IFSET(ctrl->dbglvl, DBG_REFINEINFO, rprintf(ctrl, "K-way adaptive refinement %2d [%5d, %5d], [%5d, %5d] Cut: %f, \tBalance: %6.3f\n",
        graph->level, gpwgts[idxamin(nparts, gpwgts)], gpwgts[idxamax(nparts, gpwgts)], badminpwgt, badmaxpwgt, graph->mincut,
        1.0*nparts*gpwgts[idxamax(nparts, gpwgts)]/(1.0*idxsum(nparts, gpwgts))));

  for (pass=0; pass<npasses; pass++) {
    if (mype == 0)
      RandomPermute(nparts, pperm, 0);
    MPI_Bcast((void *)pperm, nparts, IDX_DATATYPE, 0, ctrl->comm);

    oldcut = graph->mincut;

    FastRandomPermute(nvtxs, perm, 0);
    for (c=0; c<2; c++) {
      nlupd = nsupd = nmoves = nchanged = ntodo = 0;

      /* Go once through the vertices, and moves the one that are OK */
      for (ii=0; ii<nvtxs; ii++) {
        i = perm[ii];
        if (rinfo[i].ed >= rinfo[i].id) { /* Total ED is too high */
          degrees = rinfo[i].degrees;
          from = where[i];
          vwgt = graph->vwgt[i];
          vsize = graph->vsize[i];

          if (gpwgts[from]-vwgt < badminpwgt)
            continue;   /* This cannot be moved! */

          for (k=0; k<rinfo[i].ndegrees; k++) {
            other = degrees[k].edge;
            if (ProperSide(c, pperm[from], pperm[other]) && gpwgts[other]+vwgt <= badmaxpwgt)
              break;
          }

          if (k < rinfo[i].ndegrees) { /* You actually found one */
            for (j=k+1; j<rinfo[i].ndegrees; j++) {
              other = degrees[j].edge;
              if (ProperSide(c, pperm[from], pperm[other]) &&
                  ((degrees[j].ewgt > degrees[k].ewgt && gpwgts[other]+vwgt <= badmaxpwgt) ||
                  (degrees[j].ewgt == degrees[k].ewgt && gpwgts[other] < gpwgts[degrees[k].edge])))
                k = j;
            }

            other = degrees[k].edge;

            if (degrees[k].ewgt > rinfo[i].id || 
                (degrees[k].ewgt-rinfo[i].id > -vsize && gpwgts[other]+vwgt < 0.8*badminpwgt) ||
                (degrees[k].ewgt == rinfo[i].id && 
                 (other == mype || gpwgts[from] > badmaxpwgt || gpwgts[other] < badminpwgt || (from != mype && (gpwgts[from] - gpwgts[other] >= vwgt))))) {
/*
myprintf(ctrl, " Node: %3d, Gain: %3d, From: %d, To: %d, Vwgt: %3d, Fwgt: %3d, Twgt: %3d, Max/Min: %d %d\n", i, degrees[k].ewgt-rinfo[i].id, from, other, vwgt, gpwgts[from], gpwgts[other], badminpwgt, badmaxpwgt);
*/

              /* Update where, weight, and ID/ED information of the vertex you moved */
              where[i] = other;
              INC_DEC(lpwgts[other], lpwgts[from], vwgt);
              INC_DEC(gpwgts[other], gpwgts[from], vwgt);

              if (htable[i] == 0) {  /* make sure you do the update */
                htable[i] = 1;
                update[nlupd++] = i;
              }

              /* Put the vertices adjacent to i into the update array */
              for (j=xadj[i]; j<xadj[i+1]; j++) {
                k = ladjncy[j];
                if (htable[k] == 0) {
                  htable[k] = 1;
                  if (k<nvtxs)
                    update[nlupd++] = k;
                  else
                    supdate[nsupd++] = k;
                }
              }
              nmoves++;
              if (graph->pexadj[i+1]-graph->pexadj[i] > 0)
                changed[nchanged++] = i;
            }
          }
        }
        if (gpwgts[where[i]] > badmaxpwgt) 
          tocheck[ntodo++] = i;
      }

      MPI_Allreduce((void *)lpwgts, (void *)gpwgts, nparts, IDX_DATATYPE, MPI_SUM, ctrl->comm);
      IFSET(ctrl->dbglvl, DBG_REFINEINFO, rprintf(ctrl, "\tLeft todo: \t%d %d\t Balance: %5.3f\n", GlobalSESum(ctrl, ntodo), GlobalSEMax(ctrl, ntodo),
              1.0*nparts*gpwgts[idxamax(nparts, gpwgts)]/(1.0*idxsum(nparts, gpwgts))));

      /* Go through the vertices for the second time, and move negative gain for balance */
      ndone = 0;
      for (ii=0; ii<ntodo; ii++) {
        i = tocheck[ii];
        from = where[i];
        vwgt = graph->vwgt[i];
        if (gpwgts[from] > badmaxpwgt) {
          degrees = rinfo[i].degrees;

          for (k=0; k<rinfo[i].ndegrees; k++) {
            other = degrees[k].edge;
            if (ProperSide(c, pperm[from], pperm[other]) && gpwgts[other]+vwgt <= gpwgts[from])
              break;
          }

          if (k < rinfo[i].ndegrees) { /* You actually found one */
            ndone += vwgt;
            for (j=k+1; j<rinfo[i].ndegrees; j++) {
              other = degrees[j].edge;
              to = degrees[k].edge;
              if (ProperSide(c, pperm[from], pperm[other])) {
                if (gpwgts[other] < badminpwgt && (gpwgts[to] > badminpwgt || degrees[j].ewgt > degrees[k].ewgt))
                  k = j;
                else if (gpwgts[to] > badminpwgt &&
                  ((degrees[j].ewgt > degrees[k].ewgt && gpwgts[other]+vwgt <= gpwgts[from]) ||
                  (degrees[j].ewgt == degrees[k].ewgt && gpwgts[other] < gpwgts[to])))
                  k = j;
              }
            }

            other = degrees[k].edge;
/*
myprintf(ctrl, "*Node: %3d, Gain: %3d, From: %d, To: %d, Vwgt: %3d, Fwgt: %3d, Twgt: %3d, Max/Min: %d %d\n", i, degrees[k].ewgt-rinfo[i].id, from, other, vwgt, gpwgts[from], gpwgts[other], badminpwgt, badmaxpwgt);
*/
            /* Update where, weight, and ID/ED information of the vertex you moved */
            where[i] = other;
            INC_DEC(lpwgts[other], lpwgts[from], vwgt);
            INC_DEC(gpwgts[other], gpwgts[from], vwgt);

            if (htable[i] == 0) {  /* make sure you do the update */
              htable[i] = 1;
              update[nlupd++] = i;
            }

            /* Put the vertices adjacent to i into the update array */
            for (j=xadj[i]; j<xadj[i+1]; j++) {
              k = ladjncy[j];
              if (htable[k] == 0) {
                htable[k] = 1;
                if (k<nvtxs)
                  update[nlupd++] = k;
                else
                  supdate[nsupd++] = k;
              }
            }
            nmoves++;
            if (graph->pexadj[i+1]-graph->pexadj[i] > 0)
              changed[nchanged++] = i;
          }
        }
      }

      MPI_Allreduce((void *)lpwgts, (void *)gpwgts, nparts, IDX_DATATYPE, MPI_SUM, ctrl->comm);
      IFSET(ctrl->dbglvl, DBG_REFINEINFO, rprintf(ctrl, "\tActual done: \t%d %d\t Balance: %5.3f\n", GlobalSESum(ctrl, ndone), GlobalSEMax(ctrl, ndone),
              1.0*nparts*gpwgts[idxamax(nparts, gpwgts)]/(1.0*idxsum(nparts, gpwgts))));

      /* myprintf(ctrl, "nmoves: %d, nlupd: %d, nsupd: %d\n", nmoves, nlupd, nsupd); */

      /* Tell everybody interested what the new where[] info is for the interface vertices */
      CommChangedInterfaceData(ctrl, graph, nchanged, changed, where, swchanges, rwchanges, wspace->pv4); 


      IFSET(ctrl->dbglvl, DBG_RMOVEINFO, rprintf(ctrl, "\t[%d %d], [%5d, %5d],  [%d %d %d]\n", 
                pass, c, badminpwgt, badmaxpwgt, 
                GlobalSESum(ctrl, nmoves), GlobalSESum(ctrl, nsupd), GlobalSESum(ctrl, nlupd)));


      /*-------------------------------------------------------------
      / Time to communicate with processors to send the vertices
      / whose degrees need to be update.
      /-------------------------------------------------------------*/
      /* Issue the receives first */
      for (i=0; i<nnbrs; i++) {
        MPI_Irecv((void *)(rupdate+sendptr[i]), sendptr[i+1]-sendptr[i], IDX_DATATYPE,
                  peind[i], 1, ctrl->comm, ctrl->rreq+i);
      }

      /* Issue the sends next. This needs some preporcessing */
      for (i=0; i<nsupd; i++) {
        htable[supdate[i]] = 0;
        supdate[i] = graph->imap[supdate[i]];
      }
      iidxsort(nsupd, supdate);

      for (j=i=0; i<nnbrs; i++) {
        otherlastvtx = vtxdist[peind[i]+1];
        for (k=j; k<nsupd && supdate[k] < otherlastvtx; k++); 
        MPI_Isend((void *)(supdate+j), k-j, IDX_DATATYPE, peind[i], 1, ctrl->comm, ctrl->sreq+i);
        j = k;
      }

      /* OK, now get into the loop waiting for the send/recv operations to finish */
      MPI_Waitall(nnbrs, ctrl->rreq, ctrl->statuses);
      for (i=0; i<nnbrs; i++) 
        MPI_Get_count(ctrl->statuses+i, IDX_DATATYPE, nupds_pe+i);
      MPI_Waitall(nnbrs, ctrl->sreq, ctrl->statuses);


      /*-------------------------------------------------------------
      / Place the recieved to-be updated vertices into update[] 
      /-------------------------------------------------------------*/
      for (i=0; i<nnbrs; i++) {
        pe_updates = rupdate+sendptr[i];
        for (j=0; j<nupds_pe[i]; j++) {
          k = pe_updates[j];
          if (htable[k-firstvtx] == 0) {
            htable[k-firstvtx] = 1;
            update[nlupd++] = k-firstvtx;
          }
        }
      }


      /*-------------------------------------------------------------
      / Update the rinfo of the vertices in the update[] array
      /-------------------------------------------------------------*/
      for (ii=0; ii<nlupd; ii++) {
        i = update[ii];
        ASSERT(ctrl, htable[i] == 1);

        htable[i] = 0;

        me = where[i];
        myrinfo = rinfo+i;

        graph->lmincut -= myrinfo->ed;

        myrinfo->ndegrees = myrinfo->id = myrinfo->ed = 0;
        degrees = myrinfo->degrees;

        for (j=xadj[i]; j<xadj[i+1]; j++) {
          other = where[ladjncy[j]];
          if (me != other) {
            myrinfo->ed += adjwgt[j];

            for (k=0; k<myrinfo->ndegrees; k++) {
              if (degrees[k].edge == other) {
                degrees[k].ewgt += adjwgt[j];
                break;
              }
            }
            if (k == myrinfo->ndegrees) {
              degrees[k].edge = other;
              degrees[k].ewgt = adjwgt[j];
              myrinfo->ndegrees++;
            }
            ASSERT(ctrl, myrinfo->ndegrees <= xadj[i+1]-xadj[i]);
          }
          else {
            myrinfo->id += adjwgt[j];
          }
        }
        graph->lmincut += myrinfo->ed;
      }
    }

    graph->mincut = GlobalSESum(ctrl, graph->lmincut)/2;

    IFSET(ctrl->dbglvl, DBG_REFINEINFO, rprintf(ctrl, "\t[%5d, %5d], [%5d, %5d] Cut: %f\n", 
          gpwgts[idxamin(nparts, gpwgts)], gpwgts[idxamax(nparts, gpwgts)], badminpwgt, badmaxpwgt, graph->mincut));

    if (graph->mincut == oldcut)
      break;
  }

  MPI_Allreduce((void *)lpwgts, (void *)gpwgts, nparts, IDX_DATATYPE, MPI_SUM, ctrl->comm);
  IMfree(&update, &nupds_pe, &htable, &changed, &perm, &pperm, &tocheck, LTERM);

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->KWayTmr));
}
