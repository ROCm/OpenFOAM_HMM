/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * coarsen.c
 *
 * This file contains code that finds a matching and performs the coarsening
 *
 * George Irene
 */

#include "parmetis.h"


/*************************************************************************
* This function finds a HEM matching between local vertices only
**************************************************************************/
void LocalMatch_HEM(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace)
{
  int i, ii, j, k, nvtxs, cnvtxs, firstvtx, maxi; 
  idxtype *xadj, *ladjncy, *vtxdist;
  idxtype *perm, *match;
  realtype *adjwgt;

  graph->match_type = MATCH_LOCAL;

  IFSET(ctrl->dbglvl, DBG_TIME, MPI_Barrier(ctrl->comm));
  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->MatchTmr));

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  ladjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  vtxdist = graph->vtxdist;
  firstvtx = vtxdist[ctrl->mype];

  match = graph->match = idxmalloc(nvtxs+graph->nrecv, "HEM_Match: match");
  idxset(nvtxs, UNMATCHED, match);
  idxset(graph->nrecv, 0, match+nvtxs);  /* Easy way to handle remote vertices */


  /*************************************************************
   * Go now and find a local matching 
   *************************************************************/
  perm = wspace->indices;
  FastRandomPermute(nvtxs, perm, 1);
  cnvtxs = 0;
  for (ii=0; ii<nvtxs; ii++) {
    i = perm[ii];
    if (match[i] == UNMATCHED) {  /* Unmatched */
      maxi = -1;

      /* Find a heavy-edge matching, if the weight of the vertex is OK */
      if (graph->vwgt[i] < graph->maxvwgt) {
        for (j=xadj[i]; j<xadj[i+1]; j++) {
          if (match[ladjncy[j]] == UNMATCHED && graph->vwgt[ladjncy[j]] <= graph->maxvwgt &&
              (maxi == -1 || adjwgt[maxi] < adjwgt[j]))
            maxi = j;
        }
      }

      if (maxi != -1) {
        k = ladjncy[maxi];
        if (i <= k) {
          match[i] = firstvtx+k + KEEP_BIT;
          match[k] = firstvtx+i;
        }
        else {
          match[i] = firstvtx+k;
          match[k] = firstvtx+i + KEEP_BIT;
        }
      }
      else {
        match[i] = (firstvtx+i) + KEEP_BIT;
      }
      cnvtxs++;
    }
  }

  CommInterfaceData(ctrl, graph, match, wspace->indices, match+nvtxs);


#ifdef DEBUG_MATCH
  PrintVector2(ctrl, nvtxs, firstvtx, match, "Match1");
#endif


  if (ctrl->dbglvl&DBG_MATCHINFO) {
    PrintVector2(ctrl, nvtxs, firstvtx, match, "Match");
    myprintf(ctrl, "Cnvtxs: %d\n", cnvtxs);
    rprintf(ctrl, "Done with matching...\n");
  }

  IFSET(ctrl->dbglvl, DBG_TIME, MPI_Barrier(ctrl->comm));
  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->MatchTmr));


  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->ContractTmr));
  Local_CreateCoarseGraph(ctrl, graph, wspace, cnvtxs);
  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->ContractTmr));

}





/*************************************************************************
* This function creates the coarser graph
**************************************************************************/
void Local_CreateCoarseGraph(CtrlType *ctrl, GraphType *graph, WorkSpaceType *wspace, int cnvtxs)
{
  int i, j, k, l, nvtxs, nedges, firstvtx, cfirstvtx;
  int npes=ctrl->npes, mype=ctrl->mype;
  int cnedges, v, u;
  idxtype *xadj, *vwgt, *vsize, *ladjncy, *vtxdist, *where;
  idxtype *match, *cmap;
  idxtype *cxadj, *cvwgt, *cvsize, *cadjncy, *cvtxdist, *cwhere;
  realtype *adjwgt;
  realtype *cadjwgt;
  GraphType *cgraph;
  int mask=(1<<13)-1, htable[8192], htableidx[8192];

  nvtxs = graph->nvtxs;

  vtxdist = graph->vtxdist;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  vsize = graph->vsize;
  ladjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where = graph->where;

  match = graph->match;

  firstvtx = vtxdist[mype];

  cmap = graph->cmap = idxmalloc(nvtxs+graph->nrecv, "CreateCoarseGraph: cmap");

  /* Initialize the coarser graph */
  cgraph = CreateGraph();
  cgraph->nvtxs = cnvtxs;
  cgraph->maxvwgt = graph->maxvwgt;
  cgraph->level = graph->level+1;

  cgraph->finer = graph;
  graph->coarser = cgraph;


  /*************************************************************
  * Obtain the vtxdist of the coarser graph 
  **************************************************************/
  cvtxdist = cgraph->vtxdist = idxmalloc(npes+1, "CreateCoarseGraph: cvtxdist");
  cvtxdist[npes] = cnvtxs;  /* Use last position in the cvtxdist as a temp buffer */

  MPI_Allgather((void *)(cvtxdist+npes), 1, IDX_DATATYPE, (void *)cvtxdist, 1, IDX_DATATYPE, ctrl->comm);

  MAKECSR(i, npes, cvtxdist);

  cgraph->gnvtxs = cvtxdist[npes];

#ifdef DEBUG_CONTRACT
  PrintVector(ctrl, npes+1, 0, cvtxdist, "cvtxdist");
#endif


  /*************************************************************
  * Construct the cmap vector 
  **************************************************************/
  cfirstvtx = cvtxdist[mype];

  /* Create the cmap of what you know so far locally */
  cnvtxs = 0;
  for (i=0; i<nvtxs; i++) {
    if (match[i] >= KEEP_BIT) {
      k = match[i] - KEEP_BIT;
      if (k<firstvtx+i)
        continue;  /* i has been matched via the (k,i) side */

      cmap[i] = cfirstvtx + cnvtxs++;
      if (k != firstvtx+i) {
        cmap[k-firstvtx] = cmap[i];
        match[k-firstvtx] += KEEP_BIT;  /* Add the KEEP_BIT to simplify coding */
      }
    }
  }

  CommInterfaceData(ctrl, graph, cmap, wspace->indices, cmap+nvtxs);


#ifdef DEBUG_CONTRACT
  PrintVector(ctrl, nvtxs, firstvtx, cmap, "Cmap");
#endif



  /*************************************************************
  * Finally, create the coarser graph
  **************************************************************/
  /* Allocate memory for the coarser graph, and fire up coarsening */
  cxadj = cgraph->xadj = idxmalloc(cnvtxs+1, "CreateCoarserGraph: cxadj");
  cvwgt = cgraph->vwgt = idxmalloc(cnvtxs, "CreateCoarserGraph: cvwgt");
  if (vsize != NULL)
    cvsize = cgraph->vsize = idxmalloc(cnvtxs, "CreateCoarserGraph: cvsize");
  if (where != NULL)
    cwhere = cgraph->where = idxmalloc(cnvtxs, "CreateCoarserGraph: cwhere");
  cadjncy = idxmalloc(graph->nedges, "CreateCoarserGraph: cadjncy");
  cadjwgt = realmalloc(graph->nedges, "CreateCoarserGraph: cadjwgt");

  iset(8192, -1, htable);

  cxadj[0] = cnvtxs = cnedges = 0;
  for (i=0; i<nvtxs; i++) {
    v = firstvtx+i; 
    u = match[i]-KEEP_BIT;

    if (v > u) 
      continue;  /* I have already collapsed it as (u,v) */

    /* Collapse the v vertex first, which you know that is local */
    cvwgt[cnvtxs] = vwgt[i];
    if (vsize != NULL)
      cvsize[cnvtxs] = vsize[i];
    if (where != NULL)
      cwhere[cnvtxs] = where[i];
    nedges = 0;

    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = cmap[ladjncy[j]];
      if (k != cfirstvtx+cnvtxs) {  /* If this is not an internal edge */
        l = k&mask;
        if (htable[l] == -1) { /* Seeing this for first time */
          htable[l] = k;
          htableidx[l] = cnedges+nedges;
          cadjncy[cnedges+nedges] = k; 
          cadjwgt[cnedges+nedges++] = adjwgt[j];
        }
        else if (htable[l] == k) {
          cadjwgt[htableidx[l]] += adjwgt[j];
        }
        else { /* Now you have to go and do a search. Expensive case */
          for (l=0; l<nedges; l++) {
            if (cadjncy[cnedges+l] == k)
              break;
          }
          if (l < nedges) {
            cadjwgt[cnedges+l] += adjwgt[j];
          }
          else {
            cadjncy[cnedges+nedges] = k; 
            cadjwgt[cnedges+nedges++] = adjwgt[j];
          }
        }
      }
    }

    /* Collapse the u vertex next */
    if (v != u) { 
      u -= firstvtx;
      cvwgt[cnvtxs] += vwgt[u];
      if (vsize != NULL)
        cvsize[cnvtxs] += vsize[u];
      if (where != NULL && cwhere[cnvtxs] != where[u])
        myprintf(ctrl, "Something went wrong with the where local matching! %d %d\n", cwhere[cnvtxs], where[u]);

      for (j=xadj[u]; j<xadj[u+1]; j++) {
        k = cmap[ladjncy[j]];
        if (k != cfirstvtx+cnvtxs) {  /* If this is not an internal edge */
          l = k&mask;
          if (htable[l] == -1) { /* Seeing this for first time */
            htable[l] = k;
            htableidx[l] = cnedges+nedges;
            cadjncy[cnedges+nedges] = k; 
            cadjwgt[cnedges+nedges++] = adjwgt[j];
          }
          else if (htable[l] == k) {
            cadjwgt[htableidx[l]] += adjwgt[j];
          }
          else { /* Now you have to go and do a search. Expensive case */
            for (l=0; l<nedges; l++) {
              if (cadjncy[cnedges+l] == k)
                break;
            }
            if (l < nedges) {
              cadjwgt[cnedges+l] += adjwgt[j];
            }
            else {
              cadjncy[cnedges+nedges] = k; 
              cadjwgt[cnedges+nedges++] = adjwgt[j];
            }
          }
        }
      }
    }

    cnedges += nedges;
    for (j=cxadj[cnvtxs]; j<cnedges; j++)
      htable[cadjncy[j]&mask] = -1;  /* reset the htable */
    cxadj[++cnvtxs] = cnedges;
  }

  cgraph->nedges = cnedges;

  cgraph->adjncy = idxmalloc(cnedges, "CreateCoarserGraph: cadjncy");
  cgraph->adjwgt = realmalloc(cnedges, "CreateCoarserGraph: cadjwgt");
  idxcopy(cnedges, cadjncy, cgraph->adjncy);
  realcopy(cnedges, cadjwgt, cgraph->adjwgt);
  IMfree(&cadjncy, &cadjwgt, &graph->where, LTERM); /* Note that graph->where works fine even if it is NULL */
}
