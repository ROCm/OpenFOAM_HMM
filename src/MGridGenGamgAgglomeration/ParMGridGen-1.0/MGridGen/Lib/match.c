/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * match.c
 *
 * This file contains the matching routines coarsening
 *
 * George Irene
 */

#include "mgridgen.h"


/*************************************************************************
* This function finds a matching using the HEM heuristic
**************************************************************************/
void Match_RM(CtrlType *ctrl, GraphType *graph)
{
  int i, ii, j, k, nvtxs, cnvtxs, maxidx;
  idxtype *xadj, *vwgt, *adjncy;
  idxtype *match, *cmap, *perm;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;

  cmap = graph->cmap = idxsmalloc(nvtxs, -1, "graph->cmap");
  match = idxsmalloc(nvtxs, -1, "match");
  perm = idxmalloc(nvtxs, "perm");

  RandomPermute(nvtxs, perm, 1);

  cnvtxs = 0;
  for (ii=0; ii<nvtxs; ii++) {
     i = perm[ii];

     if (match[i] == UNMATCHED) {
       maxidx = i;

       /* Find a random matching, subject to maxvwgt constraints */
       for (j=xadj[i]; j<xadj[i+1]; j++) {
          k = adjncy[j];
          if (match[k] == UNMATCHED && vwgt[i]+vwgt[k] <= ctrl->maxsize) {
            maxidx = k;
            break;
          }
       }

       cmap[i] = cmap[maxidx] = cnvtxs++;
       match[i] = maxidx;
       match[maxidx] = i;
     }
  }

  CreateCoarseGraph(graph, cnvtxs, match, perm);

  IMfree(&match, &perm, LTERM);
}


/*************************************************************************
* This function finds a matching using the HEM heuristic
**************************************************************************/
void Match_HEM(CtrlType *ctrl, GraphType *graph)
{
  int i, ii, j, k, nvtxs, cnvtxs, maxidx, dim;
  idxtype *xadj, *vwgt, *adjncy;
  idxtype *match, *cmap, *perm, *tperm;
  realtype curwgt, maxwgt;
  realtype *vvol, *vsurf, *adjwgt, *adjwgtsum;

  dim       = ctrl->dim;
  nvtxs     = graph->nvtxs;
  xadj      = graph->xadj;
  vwgt      = graph->vwgt;
  vvol      = graph->vvol;
  vsurf     = graph->vsurf;
  adjncy    = graph->adjncy;
  adjwgt    = graph->adjwgt;
  adjwgtsum = graph->adjwgtsum;

  cmap = graph->cmap = idxsmalloc(nvtxs, -1, "cmap");
  match = idxsmalloc(nvtxs, -1, "match");

  perm = idxmalloc(nvtxs, "perm");
  tperm = idxmalloc(nvtxs, "tperm");

  RandomPermute(nvtxs, tperm, 1);
  BucketSortKeysInc(nvtxs, vwgt[iamax(nvtxs, vwgt)], vwgt, tperm, perm);
  /* RandomPermute(nvtxs, perm, 1);  */

  cnvtxs = 0;

  /* Compute a heavy-edge style matching giving preferance to small vertices */
  for (ii=0; ii<nvtxs; ii++) {
     i = perm[ii];

     if (match[i] == UNMATCHED) {
       maxidx = i;
       maxwgt = 0.0;

       /* Find a heavy-edge matching, subject to maxvwgt constraints */
       for (j=xadj[i]; j<xadj[i+1]; j++) {
          k = adjncy[j];
          curwgt = 1.0/ARATIO2(dim, vsurf[i]+vsurf[k]+adjwgtsum[i]+adjwgtsum[k]-
                   2.0*adjwgt[j], vvol[i]+vvol[k]);
          if (match[k] == UNMATCHED && vwgt[i]+vwgt[k] <= ctrl->maxsize &&
              curwgt > maxwgt) {
            maxwgt = curwgt;
            maxidx = k;
          }
       }

       cmap[i] = cmap[maxidx] = cnvtxs++;
       match[i] = maxidx;
       match[maxidx] = i;
     }
  }

  CreateCoarseGraph(graph, cnvtxs, match, perm);

  IMfree(&tperm, &perm, &match, LTERM);
}


/*************************************************************************
* This function finds a matching using the HEM heuristic
**************************************************************************/
void Match_HEM_Slow(CtrlType *ctrl, GraphType *graph)
{
  int i, ii, j, k, dim, nvtxs, cnvtxs, maxidx, nmatched;
  idxtype *xadj, *vwgt, *adjncy;
  idxtype *match, *cmap, *perm, *tperm;
  realtype curwgt, maxwgt;
  realtype *vvol, *vsurf, *adjwgt, *adjwgtsum;

  dim       = ctrl->dim;
  nvtxs     = graph->nvtxs;
  xadj      = graph->xadj;
  vwgt      = graph->vwgt;
  vvol      = graph->vvol;
  vsurf     = graph->vsurf;
  adjncy    = graph->adjncy;
  adjwgt    = graph->adjwgt;
  adjwgtsum = graph->adjwgtsum;

  cmap = graph->cmap = idxsmalloc(nvtxs, -1, "cmap");
  match = idxsmalloc(nvtxs, -1, "match");

  perm = idxmalloc(nvtxs, "perm");
  tperm = idxmalloc(nvtxs, "tperm");

  RandomPermute(nvtxs, tperm, 1);
  BucketSortKeysInc(nvtxs, vwgt[iamax(nvtxs, vwgt)], vwgt, tperm, perm);
  /* RandomPermute(nvtxs, perm, 1); */

  cnvtxs = 0;

  /* Compute a heavy-edge style matching giving preferance to small vertices */
  for (nmatched=0, ii=0; ii<nvtxs; ii++) {
     i = perm[ii];

     if (match[i] == UNMATCHED) {
       maxidx = i;
       maxwgt = 0.0;

       /* Find a heavy-edge matching, subject to maxvwgt constraints */
       if (nmatched < .25*nvtxs) {
         for (j=xadj[i]; j<xadj[i+1]; j++) {
            k = adjncy[j];
            if (match[k] == UNMATCHED) {
              curwgt = 1.0/ARATIO2(dim, vsurf[i]+vsurf[k]+adjwgtsum[i]+adjwgtsum[k]
                       -2.0*adjwgt[j], vvol[i]+vvol[k]);
              if (vwgt[i]+vwgt[k] <= ctrl->maxsize && curwgt > maxwgt) {
                maxwgt = curwgt;
                maxidx = adjncy[j];
              }
            }
         }
       }
       if (maxidx != i)
         nmatched++;

       cmap[i] = cmap[maxidx] = cnvtxs++;
       match[i] = maxidx;
       match[maxidx] = i;
     }
  }

  CreateCoarseGraph(graph, cnvtxs, match, perm);

  IMfree(&tperm, &perm, &match, LTERM);
}


/*************************************************************************
* This function finds a matching using the HEM heuristic
**************************************************************************/
void Match_HEM_Slow_Restricted(CtrlType *ctrl, GraphType *graph)
{
  int i, ii, j, k, dim, nvtxs, cnvtxs, maxidx, nmatched;
  idxtype *xadj, *vwgt, *adjncy, *where;
  idxtype *match, *cmap, *perm;
  realtype curwgt, maxwgt;
  realtype *vvol, *vsurf, *adjwgt, *adjwgtsum;

  dim       = ctrl->dim;
  nvtxs     = graph->nvtxs;
  xadj      = graph->xadj;
  vwgt      = graph->vwgt;
  vvol      = graph->vvol;
  vsurf     = graph->vsurf;
  adjncy    = graph->adjncy;
  adjwgt    = graph->adjwgt;
  adjwgtsum = graph->adjwgtsum;
  where     = graph->where;

  cmap = graph->cmap = idxsmalloc(nvtxs, -1, "cmap");
  match = idxsmalloc(nvtxs, -1, "match");

  perm = idxmalloc(nvtxs, "perm");
  RandomPermute(nvtxs, perm, 1);

  cnvtxs = 0;

  /* Compute a heavy-edge style matching giving preferance to small vertices */
  for (nmatched=0, ii=0; ii<nvtxs; ii++) {
     i = perm[ii];

     if (match[i] == UNMATCHED) {
       maxidx = i;
       maxwgt = 0.0;

       /* Find a heavy-edge matching, subject to maxvwgt constraints */
       if (nmatched < .3*nvtxs) {
         for (j=xadj[i]; j<xadj[i+1]; j++) {
            k = adjncy[j];
            if (where[i] != where[k])
              continue;  /* perform a restricted matching */

            curwgt = 1.0/ARATIO2(dim, vsurf[i]+vsurf[k]+adjwgtsum[i]+adjwgtsum[k]
                     -2.0*adjwgt[j], vvol[i]+vvol[k]);
            if (match[k] == UNMATCHED && vwgt[i]+vwgt[k] <= ctrl->maxsize &&
                curwgt > maxwgt) {
              maxwgt = curwgt;
              maxidx = k;
            }
         }
       }
       if (maxidx != i)
         nmatched++;

       cmap[i] = cmap[maxidx] = cnvtxs++;
       match[i] = maxidx;
       match[maxidx] = i;
     }
  }

  CreateCoarseGraph(graph, cnvtxs, match, perm);

  IMfree(&perm, &match, LTERM);
}


/*************************************************************************
* This function finds a matching using the HEM heuristic
**************************************************************************/
void Match_HEM_True(CtrlType *ctrl, GraphType *graph)
{
  int i, ii, j, k, dim, nvtxs, cnvtxs, ncand;
  idxtype *xadj, *vwgt, *adjncy;
  idxtype *match, *cmap, *perm;
  realtype *vvol, *vsurf, *adjwgt, *adjwgtsum;
  FKeyValueType *cand;

  dim       = ctrl->dim;
  nvtxs     = graph->nvtxs;
  xadj      = graph->xadj;
  vwgt      = graph->vwgt;
  vvol      = graph->vvol;
  vsurf     = graph->vsurf;
  adjncy    = graph->adjncy;
  adjwgt    = graph->adjwgt;
  adjwgtsum = graph->adjwgtsum;

  cmap = graph->cmap = idxsmalloc(nvtxs, -1, "cmap");
  match = idxsmalloc(nvtxs, -1, "match");

  perm = idxmalloc(nvtxs, "perm");
  RandomPermute(nvtxs, perm, 1);

  cand = (FKeyValueType *)IMmalloc((xadj[nvtxs]/2)*sizeof(FKeyValueType), "cand");

  /* insert the vertices according to their aspect ratio */
  for (ncand=0, ii=0; ii<nvtxs; ii++) {
     i = perm[ii];
     for (j=xadj[i]; j<xadj[i+1]; j++) {
        k = adjncy[j];
        if (k > i || vwgt[i] + vwgt[k] > ctrl->maxsize)
          continue;

        cand[ncand].val1 = i;
        cand[ncand].val2 = k;
        cand[ncand].key = ARATIO2(dim, vsurf[i]+vsurf[k]+adjwgtsum[i]+adjwgtsum[k]
                                  -2.0*adjwgt[j], vvol[i]+vvol[k]);
        ncand++;
     }
  }

  ifkeysort(ncand, cand);

  /* Compute heaviest style matching */
  idxset(nvtxs, -1, perm);
  for (cnvtxs=0, ii=0; ii<ncand; ii++) {
     if (cnvtxs > .25*nvtxs)
       break;

     i = cand[ii].val1;
     k = cand[ii].val2;

     if (match[i] == UNMATCHED && match[k] == UNMATCHED) {
       perm[cnvtxs] = i;
       perm[nvtxs-cnvtxs-1] = k;
       cmap[i] = cmap[k] = cnvtxs++;
       match[i] = k;
       match[k] = i;
     }
  }

  /* take care of the unmatched vertices */
  for (i=0; i<nvtxs; i++) {
     if (match[i] == UNMATCHED) {
       perm[cnvtxs] = i;
       cmap[i] = cnvtxs++;
       match[i] = i;
     }
  }

  CreateCoarseGraph(graph, cnvtxs, match, perm);

  IMfree(&cand, &perm, &match, LTERM);
}
