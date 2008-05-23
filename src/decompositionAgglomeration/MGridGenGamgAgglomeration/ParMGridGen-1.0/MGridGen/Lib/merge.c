/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * merge.c
 *
 * This file contains the routines for enforcing the
 * minsize and maxsize parameters
 *
 * Irene
 */

#include "mgridgen.h"

/*************************************************************************
*  This function is the entry point of merging/contributing
***************************************************************************/
void Cycle(CtrlType *ctrl, GraphType *graph, int npasses)
{
   int pass;


   for (pass = 0; pass < npasses; pass++) {
      IFSET(ctrl->dbglvl, DBG_CONTR, printf("Contribute:PASS %d\n",pass));

      Merge(ctrl, graph, npasses);
      Contribute(ctrl, graph, npasses);

      if (graph->nmoves == -1)
        break;
   }

   IMfree(&graph->pwgts, &graph->pvol, &graph->psurf, LTERM);
   ComputeKWayPartitionParams(ctrl, graph);
   Random_KWayMultiObjRefine(ctrl, graph, npasses);
}


/*************************************************************************
*  This function is the entry point of merging
***************************************************************************/
void Merge(CtrlType *ctrl, GraphType *graph, int npasses)
{
   int pass;

   npasses = 2;
   for (pass=0; pass < npasses; pass++) {
      IFSET(ctrl->dbglvl, DBG_MERGE, printf("Merge: Pass %d\n",pass)); 

      switch (ctrl->RType) {
        case REFINE_AR:
          Merge_ARatio(ctrl, graph);
          break;
        case REFINE_WAR:
          Merge_WeightARatio(ctrl, graph); 
          break;
        case REFINE_SCUT:
        case REFINE_MINMAXAVAR:
          Merge_MultiObj(ctrl, graph);
          break;
        case REFINE_MINMAXAR:
          Merge_MinMaxARatio(ctrl, graph);
          break;
        case REFINE_MULTIOBJECTIVE:
        case REFINE_MULTIOBJECTIVE2:
          Merge_MultiObj(ctrl, graph);
          break;
        default:
          errexit("Unknown RType of %d\n", ctrl->RType);
      }

      if (graph->nmoves == 0)
         break;
   }
}


/*************************************************************************
*  This function is the entry point of contributing
***************************************************************************/
void Contribute(CtrlType *ctrl, GraphType *graph, int npasses)
{
   int pass;

   for (pass = 0; pass < npasses; pass++) {
      IFSET(ctrl->dbglvl, DBG_CONTR, printf("Contribute:PASS %d\n",pass)); 

      switch (ctrl->RType) {
        case REFINE_AR:
	  Contribute_ARatio(ctrl, graph);
          break;
        case REFINE_WAR:
	  Contribute_WeightARatio(ctrl, graph);
	  break;
        case REFINE_SCUT:
        case REFINE_MINMAXAVAR:
          Contribute_MultiObj(ctrl, graph);
          break;
        case REFINE_MINMAXAR:
          Contribute_MinMaxARatio(ctrl, graph);
          break;
        case REFINE_MULTIOBJECTIVE:
        case REFINE_MULTIOBJECTIVE2:
          Contribute_MultiObj(ctrl, graph);
          break;
        default:
          errexit("Unknown RType of %d\n", ctrl->RType);
      }

      if (graph->nmoves == -1)
        break;
   }
}


/*****************************************************************************
* This function merges small fused elements.
* The objective is to directly minimize the sum of the aspect ratios
******************************************************************************/
void Merge_ARatio(CtrlType *ctrl, GraphType *graph)
{
  int i, j, l, last, jbest, me, to, el, nsels, dim, nmoves, ndegrees;
  int nvtxs, nparts, minsize, maxsize;
  idxtype *xadj, *vwgt, *adjncy, *where, *pwgts;
  idxtype *phtable, *ptarget;
  idxtype *ind, *ptr;
  realtype old, new, best, OldToAR, OldFromAR,  NewToAR;
  realtype *vvol, *vsurf, *adjwgt, *pvol, *psurf;
  realtype *degrees;
  idxKeyValueType *elpairs;

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  vwgt   = graph->vwgt;
  vvol   = graph->vvol;
  vsurf  = graph->vsurf;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where  = graph->where;

  dim     = ctrl->dim;
  minsize = ctrl->minsize;
  maxsize = ctrl->maxsize;
  nparts  = ctrl->nparts;

  /* Create a part-to-vertex mapping; compute pwgts, psurf, pvol */
  ptr   = idxsmalloc(nparts+1, 0, "Merge_ARatio: ptr");
  ind   = idxmalloc(nvtxs, "Merge_ARatio: ind");
  pwgts = idxsmalloc(nparts, 0, "Merge_ARatio: pwgts");
  pvol  = realsmalloc(nparts, 0.0, "Merge_ARatio: pvol");
  psurf = realsmalloc(nparts, 0.0, "Merge_ARatio: psurf");

  for (i = 0; i < nvtxs; i++) {
     me = where[i];
     ptr[me]++;
     pwgts[me] += vwgt[i];
     pvol[me] += vvol[i];
     psurf[me] += vsurf[i];
     for (j = xadj[i]; j < xadj[i+1]; j++)
        if (where[adjncy[j]] != me)
          psurf[me] += adjwgt[j];
  }

  MAKECSR(i, nparts, ptr);
  for (i =0; i<nvtxs; i++) {
     me = where[i];
     ind[ptr[me]] = i;
     ptr[me]++;
  }

  for (i = nparts; i > 0; i--)
     ptr[i] = ptr[i-1];
  ptr[0] = 0;

  /* Find elements that are < minsize and sort them according to pwgts */
  /* nsels = number of small partitions */
  elpairs =(idxKeyValueType *) IMmalloc(nvtxs*sizeof(idxKeyValueType), "Merge_ARatio: elpairs");
  for (nsels=i=0; i < nparts; i++) 
     if (pwgts[i] < minsize) {
       elpairs[nsels].key = pwgts[i];
       elpairs[nsels++].val = i;
     }
  idxkeysort(nsels, elpairs);

  IFSET(ctrl->dbglvl, DBG_MERGE,
        printf("===== nsels = %d ===== nparts = %d =====\n",nsels, nparts));

  ptarget = idxmalloc(nparts, "Merge_ARatio: ptarget");
  degrees = realsmalloc(nparts, 0, "Merge_ARatio: degrees");
  phtable = idxsmalloc(nparts, -1, "Merge_ARatio: phtable");

  /* Determine the connectivity of the undersized fused elements */
  for (nmoves = 0, el = nsels-1; el >= 0; el--) {
     for (ndegrees = 0, me = elpairs[el].val, l = ptr[me]; l < ptr[me+1]; l++) {
        for (i = ind[l], j = xadj[i]; j < xadj[i+1]; j++) {
           to = where[adjncy[j]];

           if (to == me || pwgts[to]+pwgts[me] > maxsize || pwgts[to]==0) 
             continue;

           if (phtable[to] == -1) {        /* connection is not created yet */
             ptarget[ndegrees] = to;
             degrees[ndegrees] = adjwgt[j];
             phtable[to] = ndegrees++;
           }
           else                            /* connection is already there */
             degrees[phtable[to]] += adjwgt[j];
        }
     }

     /* Determine which of the ndegrees moves is the best */
     if (ndegrees > 0) {
       j = 0;
       to = ptarget[j];

       OldFromAR = ARATIO(dim, psurf[me], pvol[me]);
       OldToAR   = ARATIO(dim, psurf[to], pvol[to]);
       NewToAR   = ARATIO(dim, psurf[me] + psurf[to] - 2.0*degrees[j], pvol[me] + pvol[to]);

       old = OldFromAR + OldToAR;
       new = NewToAR;
 
       jbest = j;
       best = old-new;

       for (j = 1; j < ndegrees; j++) {
          to = ptarget[j];

	  OldToAR = ARATIO(dim, psurf[to], pvol[to]);
          NewToAR = ARATIO(dim, psurf[me] + psurf[to] - 2.0*degrees[j], pvol[me] + pvol[to]);

          old = OldFromAR + OldToAR;
          new = NewToAR;

          if (best < old-new) {
             jbest = j;
             best = old-new;
          }
       }

       /* Merge now */
       nmoves++;
       to = ptarget[jbest];
       IFSET(ctrl->dbglvl, DBG_MERGE,
             printf("Ndeg=%d Merge %d and %d best=%f\n", ndegrees, me, to, best));

       for (l = ptr[me]; l < ptr[me+1]; l++) {
          i = ind[l];
          where[i] = to;
       }
       INC_DEC(pwgts[to], pwgts[me], pwgts[me]);
       INC_DEC(pvol[to], pvol[me], pvol[me]);
       psurf[to] = psurf[to] + psurf[me] - 2.0*degrees[jbest];
       psurf[me] = 0;
     }

     for (j = 0; j < ndegrees; j++)
        phtable[ptarget[j]] = -1;
  }

  /* Remove empty partitions */
  for (i = 0; i < nvtxs; i++) {
     elpairs[i].key = where[i];
     elpairs[i].val = i;
  }
  idxkeysort(nvtxs, elpairs);

  for (last = 0, where[elpairs[0].val] = last, i = 1; i < nvtxs; i++) {
     if (elpairs[i].key > elpairs[i-1].key)
       last++;
     where[elpairs[i].val] = last;
  }
  ctrl->nparts  = last+1;

  graph->nmoves = nmoves;

  IMfree(&ptr, &ind, &psurf, &pvol, &pwgts, &elpairs, &ptarget, &degrees, &phtable, LTERM);
}


/*****************************************************************************
* This function merges small fused elements.
* The objective is to directly minimize the weighted sum of the aspect ratios
******************************************************************************/
void Merge_WeightARatio(CtrlType *ctrl, GraphType *graph)
{
  int i, j, l;
  int last, jbest, me, to, nparts, el, nsels, dim, nmoves, ndegrees;
  int nvtxs, minsize, maxsize;
  idxtype *xadj, *vwgt, *adjncy, *where, *pwgts;
  idxtype *phtable, *ptarget;
  idxtype *ind, *ptr;
  realtype old, new, best, OldToAR, OldFromAR,  NewToAR;
  realtype *vvol, *vsurf, *adjwgt, *pvol, *psurf;
  realtype *degrees;
  idxKeyValueType *elpairs;


  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  vwgt   = graph->vwgt;
  vvol   = graph->vvol;
  vsurf  = graph->vsurf;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where  = graph->where;

  dim     = ctrl->dim;
  minsize = ctrl->minsize;
  maxsize = ctrl->maxsize;
  nparts  = ctrl->nparts;

  /* Create a part-to-vertex mapping; compute pwgts, psurf, pvol */
  ptr   = idxsmalloc(nparts+1, 0, "Merge_WeightARatio: ptr");
  ind   = idxmalloc(nvtxs, "Merge_WeightARatio: ind");
  pwgts = idxsmalloc(nparts, 0, "Merge_WeightARatio: pwgts");
  pvol  = realsmalloc(nparts, 0.0, "Merge_WeightARatio: pvol");
  psurf = realsmalloc(nparts, 0.0, "Merge_WeightARatio: psurf");

  for (i = 0; i < nvtxs; i++) {
     me = where[i];
     ptr[me]++;
     pwgts[me] += vwgt[i];
     pvol[me] += vvol[i];
     psurf[me] += vsurf[i];
     for (j = xadj[i]; j < xadj[i+1]; j++)
        if (where[adjncy[j]] != me)
          psurf[me] += adjwgt[j];
  }

  MAKECSR(i, nparts, ptr);
  for (i =0; i<nvtxs; i++) {
     me = where[i];
     ind[ptr[me]] = i;
     ptr[me]++;
  }

  for (i = nparts; i > 0; i--)
     ptr[i] = ptr[i-1];
  ptr[0] = 0;

  /* Find elements that are < minsize and sort them according to pwgts */
  /* nsels = number of small partitions */
  elpairs =(idxKeyValueType *) IMmalloc(nvtxs*sizeof(idxKeyValueType), "Merge_WeightARatio: elpairs");
  for (nsels=i=0; i < nparts; i++) 
     if (pwgts[i] < minsize) {
       elpairs[nsels].key = pwgts[i];
       elpairs[nsels++].val = i;
     }
  idxkeysort(nsels, elpairs);

  IFSET(ctrl->dbglvl, DBG_MERGE,
        printf("===== nsels = %d ===== nparts = %d =====\n",nsels, nparts));

  ptarget = idxmalloc(nparts, "Merge_WeightARatio: ptarget");
  degrees = realsmalloc(nparts, 0, "Merge_WeightARatio: degrees");
  phtable = idxsmalloc(nparts, -1, "Merge_WeightARatio: phtable");

  /* Determine the connectivity of the undersized fused elements */
  for (nmoves = 0, el = nsels-1; el >= 0; el--) {
     for (ndegrees = 0, me = elpairs[el].val, l = ptr[me]; l < ptr[me+1]; l++) {
        for (i = ind[l], j = xadj[i]; j < xadj[i+1]; j++) {
           to = where[adjncy[j]];

           if (to == me || pwgts[to]+pwgts[me] > maxsize || pwgts[to]==0) 
             continue;

           if (phtable[to] == -1) {        /* connection is not created yet */
             ptarget[ndegrees] = to;
             degrees[ndegrees] = adjwgt[j];
             phtable[to] = ndegrees++;
           }
           else                            /* connection is already there */
             degrees[phtable[to]] += adjwgt[j];
        }
     }

     /* Determine which of the ndegrees moves is the best */
     if (ndegrees > 0) {
       j = 0;
       to = ptarget[j];

       OldFromAR = ARATIO(dim, psurf[me], pvol[me]) * pwgts[me];
       OldToAR   = ARATIO(dim, psurf[to], pvol[to]);
       NewToAR   = ARATIO(dim, psurf[me] + psurf[to] -2.0*degrees[j], pvol[me] + pvol[to]);

       old = OldFromAR + OldToAR*pwgts[to];
       new = NewToAR*(pwgts[me]+pwgts[to]);
 
       jbest = j;
       best = old-new;

       for (j = 1; j < ndegrees; j++) {
          to = ptarget[j];

	  OldToAR = ARATIO(dim, psurf[to], pvol[to]);
          NewToAR = ARATIO(dim, psurf[me] + psurf[to] -2.0*degrees[j], pvol[me] + pvol[to]);

          old = OldFromAR + OldToAR*pwgts[to];
          new = NewToAR*(pwgts[me]+pwgts[to]);

          if (best < old-new) {
             jbest = j;
             best = old-new;
          }
       }

       /* Merge now */
       nmoves++;
       to = ptarget[jbest];
       IFSET(ctrl->dbglvl, DBG_MERGE,
             printf("Ndeg=%d Merge %d and %d best=%f\n", ndegrees, me, to, best));

       for (l = ptr[me]; l < ptr[me+1]; l++) {
          i = ind[l];
          where[i] = to;
       }
       INC_DEC(pwgts[to], pwgts[me], pwgts[me]);
       INC_DEC(pvol[to], pvol[me], pvol[me]);
       psurf[to] = psurf[to] + psurf[me] - 2.0*degrees[jbest];
       psurf[me] = 0;
     }

     for (j = 0; j < ndegrees; j++)
        phtable[ptarget[j]] = -1;
  }

  /* Remove empty partitions */
  for (i = 0; i < nvtxs; i++) {
     elpairs[i].key = where[i];
     elpairs[i].val = i;
  }
  idxkeysort(nvtxs, elpairs);

  for (last = 0, where[elpairs[0].val] = last, i = 1; i < nvtxs; i++) {
     if (elpairs[i].key > elpairs[i-1].key)
       last++;
     where[elpairs[i].val] = last;
  }
  ctrl->nparts  = last+1;

  graph->nmoves = nmoves;

  IMfree(&ptr, &ind, &psurf, &pvol, &pwgts, &elpairs, &ptarget, &degrees, &phtable, LTERM);
}


/*************************************************************************
* This function merges small fused elements 
* The objective is to directly minimize the maximum aspect ratio.
**************************************************************************/
void Merge_MinMaxARatio(CtrlType *ctrl, GraphType *graph)
{
  int i, j, l, el;
  int last, me, to, nsels, ndegrees, nmoves, pmax;
  int jbest; 
  int dim, nvtxs, nparts, minsize, maxsize;
  idxtype *xadj, *vwgt, *adjncy, *where, *pwgts;
  idxtype *phtable, *ptarget;
  idxtype *ind, *ptr;
  realtype best, new, maxar, NewToAR;
  realtype *vvol, *vsurf, *adjwgt, *pvol, *psurf;
  realtype *degrees;
  idxKeyValueType *elpairs;


  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  vwgt   = graph->vwgt;
  vvol   = graph->vvol;
  vsurf  = graph->vsurf;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where  = graph->where;

  dim     = ctrl->dim;
  minsize = ctrl->minsize;
  maxsize = ctrl->maxsize;
  nparts  = ctrl->nparts;

  /* Create a part-to-vertex mapping; compute pwgts, psurf, pvol */
  ptr   = idxsmalloc(nparts+1, 0, "Merge_MinMaxARatio: ptr");
  ind   = idxmalloc(nvtxs, "Merge_MinMaxARatio: ind");
  pwgts = idxsmalloc(nparts, 0, "Merge_MinMaxARatio: pwgts");
  pvol  = realsmalloc(nparts, 0.0, "Merge_MinMaxARatio: pvol");
  psurf = realsmalloc(nparts, 0.0, "Merge_MinMaxARatio: psurf");

  for (i = 0; i < nvtxs; i++) {
     me = where[i];
     ptr[me]++;
     pwgts[me] += vwgt[i];
     pvol[me] += vvol[i];
     psurf[me] += vsurf[i];
     for (j = xadj[i]; j < xadj[i+1]; j++)
        if (where[adjncy[j]] != me)
          psurf[me] += adjwgt[j];
  }

  MAKECSR(i, nparts, ptr);
  for (i = 0; i < nvtxs; i++) {
     me = where[i];
     ind[ptr[me]] = i;
     ptr[me]++;
  }

  for (i = nparts; i>0; i--)
     ptr[i] = ptr[i-1];
  ptr[0] = 0;

  /* Find elements that are < minsize and sort them according to pwgts */
  /* nsels = number of small elements */
  elpairs =(idxKeyValueType *) IMmalloc(nvtxs*sizeof(idxKeyValueType), "Merge_MinMaxARatio: elpairs");
  for (nsels=i=0; i < nparts; i++)
     if (pwgts[i] < minsize) {
       elpairs[nsels].key = pwgts[i];
       elpairs[nsels++].val = i;
     }
  idxkeysort(nsels, elpairs);

  IFSET(ctrl->dbglvl, DBG_MERGE,
        printf("===== nsels = %d ===== nparts = %d =====\n", nsels, nparts));

  ptarget = idxmalloc(nparts, "Merge_MinMaxARatio: ptarget");
  degrees = realsmalloc(nparts, 0, "Merge_MinMaxARatio: degrees");
  phtable = idxsmalloc(nparts, -1, "Merge_MinMaxARatio: phtable"); 

  /* Determine the domain that has the maximum aspect ratio */
  pmax = 0;
  maxar = ARATIO(dim, psurf[0], pvol[0]);
  for (i = 1; i < nparts; i++) {
     new = ARATIO(dim, psurf[i], pvol[i]);
     if (new > maxar) {
       maxar = new;
       pmax = i;
     }
  }

  /* Determine the connectivity of the undersized fused elements */
  for (nmoves = 0, el = nsels-1; el >= 0; el--) {
     for (ndegrees = 0, me = elpairs[el].val, l = ptr[me]; l < ptr[me+1]; l++) {
        for (i = ind[l], j = xadj[i]; j < xadj[i+1]; j++) {
           to = where[adjncy[j]];
       
           if (to == me || pwgts[to]+pwgts[me] > maxsize || pwgts[to]==0) 
             continue;

           if (phtable[to] == -1) {        /* connection is not created yet */
             ptarget[ndegrees] = to;
             degrees[ndegrees] = adjwgt[j]; 
             phtable[to] = ndegrees++;
           }
           else                            /* connection is already there */
             degrees[phtable[to]] += adjwgt[j];
        }
     }

     /* Determine which of the ndegrees moves is the best */
     if (ndegrees > 0) {
       j = 0;
       to = ptarget[j];

       NewToAR   = ARATIO(dim, psurf[me]+psurf[to]-2.0*degrees[j], pvol[me]+pvol[to]);

       jbest = j;
       best = NewToAR;

       for (j = 1; j < ndegrees; j++) {
          to = ptarget[j];

          NewToAR   = ARATIO(dim, psurf[me]+psurf[to]-2.0*degrees[j], pvol[me]+pvol[to]);

          /* Check if it increases the aspect ratio */
          if (NewToAR < best) {
            jbest = j;
            best = NewToAR;
          }
       }

       /* Merge now */
       nmoves++;
       to = ptarget[jbest];

       for (l = ptr[me]; l < ptr[me+1]; l++) {
          i = ind[l];
          where[i] = to;
       }
       INC_DEC(pwgts[to], pwgts[me], pwgts[me]);
       INC_DEC(pvol[to], pvol[me], pvol[me]);
       psurf[to] = psurf[to] + psurf[me] - 2.0*degrees[jbest];
       psurf[me] = 0;

       /* find the new maximum aspect ratio */
       if (best > maxar || to == pmax || me == pmax )
         for (pmax=-1, i = 0; i < nparts; i++)
            if (psurf[i] != 0) {
              new = ARATIO(dim, psurf[i], pvol[i]);
              if (pmax == -1 || new > maxar) {
                maxar = new;
                pmax = i;
              }
            }

       IFSET(ctrl->dbglvl, DBG_MERGE,
             printf("Ndeg=%d Merge %d with %d jbest=%d best=%f pmax=%d maxar=%f\n",
                    ndegrees, me, to, jbest, best, pmax, maxar));

     }

     for (j = 0; j < ndegrees; j++)
        phtable[ptarget[j]] = -1;
  }

  /* Remove empty elements */
  for (i = 0; i < nvtxs; i++) {
     elpairs[i].key = where[i];
     elpairs[i].val = i;
  }
  idxkeysort(nvtxs, elpairs);

  for (last = 0, where[elpairs[0].val] = last, i = 1; i < nvtxs; i++) {
     if (elpairs[i].key > elpairs[i-1].key)
       last++;
     where[elpairs[i].val] = last;
  }
  ctrl->nparts  = last+1;

  graph->nmoves = nmoves;

  IMfree(&ptr, &ind, &psurf, &pvol, &pwgts, &elpairs, &ptarget, &degrees, &phtable, LTERM);
}


/*************************************************************************
* This function merges small fused elements.
* The objective is to directly minimize the multiple objective aspect ratio
* In this case RType = REFINE_MINMAXAR and RType =  REFINE_WAR
**************************************************************************/
void Merge_MultiObj(CtrlType *ctrl, GraphType *graph)
{
  int i, j, l, el;
  int last, me, to, nsels, ndegrees, nmoves, pmax;
  int jbest, jbest1, jbest2;
  int dim, nvtxs, nparts, minsize, maxsize;
  idxtype *xadj, *vwgt, *adjncy, *where, *pwgts;
  idxtype *phtable, *ptarget;
  idxtype *ind, *ptr;
  realtype best, best1, best2, old, new, maxar, OldToAR, OldFromAR,  NewToAR;
  realtype *vvol, *vsurf, *adjwgt, *pvol, *psurf;
  realtype *degrees;
  idxKeyValueType *elpairs;


  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  vwgt   = graph->vwgt;
  vvol   = graph->vvol;
  vsurf  = graph->vsurf;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where  = graph->where;

  dim     = ctrl->dim;
  minsize = ctrl->minsize;
  maxsize = ctrl->maxsize;
  nparts  = ctrl->nparts;

  /* Create a part-to-vertex mapping; compute pwgts, psurf, pvol */
  ptr   = idxsmalloc(nparts+1, 0, "Merge_MultiObj: ptr");
  ind   = idxmalloc(nvtxs, "Merge_MultiObj: ind");
  pwgts = idxsmalloc(nparts, 0, "Merge_MultiObj: pwgts");
  pvol  = realsmalloc(nparts, 0.0, "Merge_MultiObj: pvol");
  psurf = realsmalloc(nparts, 0.0, "Merge_MultiObj: psurf");

  for (i = 0; i < nvtxs; i++) {
     me = where[i];
     ptr[me]++;
     pwgts[me] += vwgt[i];
     pvol[me] += vvol[i];
     psurf[me] += vsurf[i];
     for (j = xadj[i]; j < xadj[i+1]; j++)
        if (where[adjncy[j]] != me)
          psurf[me] += adjwgt[j];
  }

  MAKECSR(i, nparts, ptr);
  for (i = 0; i < nvtxs; i++) {
     me = where[i];
     ind[ptr[me]] = i;
     ptr[me]++;
  }

  for (i = nparts; i > 0; i--)
     ptr[i] = ptr[i-1];
  ptr[0] = 0;

  /* Find elements that are < minsize and sort them according to pwgts */
  /* nsels = number of small elements */
  elpairs =(idxKeyValueType *) IMmalloc(nvtxs*sizeof(idxKeyValueType), "Merge_MultiObj: elpairs");
  for (nsels=i=0; i < nparts; i++) 
     if (pwgts[i] < minsize) {
       elpairs[nsels].key = pwgts[i];
       elpairs[nsels++].val = i;
     }
  idxkeysort(nsels, elpairs);

  IFSET(ctrl->dbglvl, DBG_MERGE,
        printf("===== nsels = %d ===== nparts = %d =====\n",nsels, nparts));

  ptarget = idxmalloc(nparts, "Merge_MultiObj: ptarget");
  degrees = realsmalloc(nparts, 0, "Merge_MultiObj: degrees");
  phtable = idxsmalloc(nparts, -1, "Merge_MultiObj: phtable"); 

  /* Determine the domain that has the maximum aspect ratio */
  pmax = 0;
  maxar = ARATIO(dim, psurf[0], pvol[0]);
  for (i = 1; i < nparts; i++) {
     new = ARATIO(dim, psurf[i], pvol[i]);
     if (new > maxar) {
       maxar = new;
       pmax = i;
     }
  }

  /* Determine the connectivity of the undersized fused elements */
  for (nmoves=0, el = nsels-1; el >= 0; el--) {
     for (ndegrees = 0, me = elpairs[el].val, l = ptr[me]; l < ptr[me+1]; l++) {
        for (i = ind[l], j = xadj[i]; j < xadj[i+1]; j++) {
           to = where[adjncy[j]];
       
           if (to == me || pwgts[to]+pwgts[me] > maxsize || pwgts[to]==0) 
             continue;

           if (phtable[to] == -1) {        /* connection is not created yet */
             ptarget[ndegrees] = to;
             degrees[ndegrees] = adjwgt[j]; 
             phtable[to] = ndegrees++;
           }
           else                            /* connection is already there */
             degrees[phtable[to]] += adjwgt[j];
        }
     }

     /* Determine which of the ndegrees moves is the best */
     if (ndegrees > 0) {
       j = 0;
       to = ptarget[j];

       OldFromAR = ARATIO(dim, psurf[me], pvol[me]) * pwgts[me]; 
       NewToAR   = ARATIO(dim, psurf[me]+psurf[to]-2.0*degrees[j], pvol[me]+pvol[to]);

       jbest1 = j;
       best1 = NewToAR;

       if (NewToAR <= maxar) {
         OldToAR   = ARATIO(dim, psurf[to], pvol[to]);

         old = OldFromAR + OldToAR*pwgts[to];
         new = NewToAR*(pwgts[me]+pwgts[to]);

         jbest2 = jbest1;
         best2 =  old - new;
       }
       else 
         jbest2 = -1;

       for (j = 1; j < ndegrees; j++) {
          to = ptarget[j];

          NewToAR   = ARATIO(dim, psurf[me]+psurf[to]-2.0*degrees[j], pvol[me]+pvol[to]);

          /* Check if it increases the aspect ratio */
          if (NewToAR < best1) {
            jbest1 = j;
            best1 = NewToAR;
          }

          /* If first objective is OK, check second one */
          if (NewToAR <= maxar) {
            OldToAR   = ARATIO(dim, psurf[to], pvol[to]);

            old = OldFromAR + OldToAR*pwgts[to];
            new = NewToAR*(pwgts[me]+pwgts[to]);
 
            if ( (jbest2 == -1) || (best2 < old-new) ) {
              jbest2 = j;
              best2 =  old - new;
            }
          }

       }

       if (jbest2 != -1) {
         jbest = jbest2;
         best = best2;
       }
       else {
         jbest = jbest1;
         best = best1;
       }

       /* Merge now */
       nmoves++;
       to = ptarget[jbest];

       for (l = ptr[me]; l < ptr[me+1]; l++) {
          i = ind[l];
          where[i] = to;
       }
       INC_DEC(pwgts[to], pwgts[me], pwgts[me]);
       INC_DEC(pvol[to], pvol[me], pvol[me]);
       psurf[to] = psurf[to] + psurf[me] - 2.0*degrees[jbest];
       psurf[me] = 0;

       /* find the new maximum aspect ratio */
       if (jbest2 == -1 || to == pmax || me == pmax )
         for (pmax=-1, i = 0; i < nparts; i++)
            if (psurf[i] != 0) {
              new = ARATIO(dim, psurf[i], pvol[i]);
              if (pmax == -1 || new > maxar) {
                maxar = new;
                pmax = i;
              }
            }

       IFSET(ctrl->dbglvl, DBG_MERGE,
             printf("Ndeg=%d Merge %d with %d jbest2=%d best=%f pmax=%d maxar=%f\n",
                    ndegrees, me, to, jbest2, best, pmax, maxar));
     }

     for (j = 0; j < ndegrees; j++)
        phtable[ptarget[j]] = -1;
  }

  /* Remove empty elements */
  for (i = 0; i < nvtxs; i++) {
     elpairs[i].key = where[i];
     elpairs[i].val = i;
  }
  idxkeysort(nvtxs, elpairs);

  for (last = 0, where[elpairs[0].val] = last, i = 1; i < nvtxs; i++) {
     if (elpairs[i].key > elpairs[i-1].key)
       last++;
     where[elpairs[i].val] = last;
  }
  ctrl->nparts  = last+1;
  graph->nmoves = nmoves;

  IMfree(&ptr, &ind, &psurf, &pvol, &pwgts, &elpairs, &ptarget, &degrees, &phtable, LTERM);
}


/*************************************************************************
* This function merges small fused elements 
* TWO OBJECTIVES : MAX AR & WEIGHTED SUM (DO NOT VIOLATE MAX AR)
**************************************************************************/
void Merge0(CtrlType *ctrl, GraphType *graph)
{
  int i, j, l, last, me, to, el, nsels, ndegrees, nmoves, pmax;
  int jbest1, jbest2, jbest;
  int dim, nvtxs, nparts, minsize, maxsize;
  idxtype *xadj, *vwgt, *adjncy, *where, *pwgts;
  idxtype *phtable, *ptarget;
  idxtype *ind, *ptr;
  realtype old, new, maxar, OldFromAR, OldToAR, NewToAR;
  realtype best1, best2, best;
  realtype *vvol, *vsurf, *adjwgt, *pvol, *psurf, *degrees;
  idxKeyValueType *elpairs;


  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  vwgt   = graph->vwgt;
  vvol   = graph->vvol;
  vsurf  = graph->vsurf;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where  = graph->where;

  dim     = ctrl->dim;
  minsize = ctrl->minsize;
  maxsize = ctrl->maxsize;
  nparts  = ctrl->nparts;

  /* Create a part-to-vertex mapping; compute pwgts, psurf, pvol */
  ptr   = idxsmalloc(nparts+1, 0, "Merge: ptr");
  ind   = idxmalloc(nvtxs, "Merge: ind");
  pwgts = idxsmalloc(nparts, 0, "Merge: pwgts");
  pvol  = realsmalloc(nparts, 0.0, "Merge: pvol");
  psurf = realsmalloc(nparts, 0.0, "Merge: psurf");

  for (i = 0; i < nvtxs; i++) {
     me = where[i];
     ptr[me]++;
     pwgts[me] += vwgt[i];
     pvol[me] += vvol[i];
     psurf[me] += vsurf[i];
     for (j = xadj[i]; j < xadj[i+1]; j++)
        if (where[adjncy[j]] != me)
          psurf[me] += adjwgt[j];
  }

  MAKECSR(i, nparts, ptr);
  for (i = 0; i < nvtxs; i++) {
     me = where[i];
     ind[ptr[me]] = i;
     ptr[me]++;
  }

  for (i = nparts; i > 0; i--)
     ptr[i] = ptr[i-1];
  ptr[0] = 0;

  /* Find elements that are < minsize and sort them according to pwgts */
  /* nsels = number of small elements */
  elpairs =(idxKeyValueType *) IMmalloc(nvtxs*sizeof(idxKeyValueType), "Merge:elpairs");
  for (nsels=i=0; i < nparts; i++)
     if (pwgts[i] < minsize) {
       elpairs[nsels].key = pwgts[i];
       elpairs[nsels++].val = i;
     }
  idxkeysort(nsels, elpairs);

  IFSET(ctrl->dbglvl, DBG_MERGE,
        printf("===== nsels = %d ===== nparts = %d =====\n",nsels, nparts));

  ptarget = idxmalloc(nparts, "FusedElementGraph: fadjncy");
  degrees = realsmalloc(nparts, 0, "FusedElementGraph: fadjwgt");
  phtable = idxsmalloc(nparts, -1, "FusedElementGraph: map");

  /* Determine the domain that has the maximum aspect ratio */
  pmax = 0;
  maxar = ARATIO(dim, psurf[0], pvol[0]);
  for (i = 1; i < nparts; i++) {
     new = ARATIO(dim, psurf[i], pvol[i]);
     if (new > maxar) {
       maxar = new;
       pmax = i;
     }
  }

  /* Determine the connectivity of the undersized fused elements */
  for (nmoves = 0, el = nsels-1; el >= 0; el--) {
     for (ndegrees = 0, me = elpairs[el].val, l = ptr[me]; l < ptr[me+1]; l++) {
        for (i = ind[l], j = xadj[i]; j < xadj[i+1]; j++) {
           to = where[adjncy[j]];

           if (to == me || pwgts[to]+pwgts[me] > maxsize || pwgts[to]==0)
             continue;

             if (phtable[to] == -1) {        /* connection is not created yet */
               ptarget[ndegrees] = to;
               degrees[ndegrees] = adjwgt[j];
               phtable[to] = ndegrees++;
             }
             else                            /* connection is already there */
               degrees[phtable[to]] += adjwgt[j];
        }
     }

     /* Determine which of the ndegrees moves is the best */
     if (ndegrees > 0) {
       for (best1=best2=0.0, jbest1=jbest2=-1, j = 0; j < ndegrees; j++) {
          to = ptarget[j];

          OldFromAR = ARATIO(dim, psurf[me], pvol[me]);
          OldToAR   = ARATIO(dim, psurf[to], pvol[to]);
          NewToAR   = ARATIO(dim, psurf[me]+psurf[to]-2.0*degrees[j], pvol[me]+pvol[to]);

          /* Check first objective min(max) */

          /* Check if it increases the max aspect ratio */
          if (NewToAR > maxar)
            continue;

          /* If not... */
          /* If move involves partition with max asp ratio, do the move now */
          if (to == pmax || me == pmax) {
            jbest1 = j;
            break;
          }

          /* Else if partition with max asp ratio is not involved, do the move
             that gives best local gain */
          else {
            old = amax(OldFromAR, OldToAR);
            new = NewToAR;
            if (old-new > best1) {
              best1  = old-new;
              jbest1 = j;
            }
          }

          /* Check second objective weighted sum */
          old = OldFromAR*pwgts[me] + OldToAR*pwgts[to];
          new = NewToAR*(pwgts[me]+pwgts[to]);

          if (best2 < old-new) {
            best2 = old-new;
            jbest2 = j;
          }
       }

       if (jbest1 != -1) {
         jbest = jbest1;
         best = best1;
       }
       else if (jbest2 != -1) {
         jbest = jbest2;
         best = best2;
       }
       else
         jbest = -1;


       /* Merge now */
       if (jbest != -1) {
         nmoves++;
         to = ptarget[jbest];

         for (l = ptr[me]; l < ptr[me+1]; l++) {
            i = ind[l];
            where[i] = to;
         }
         INC_DEC(pwgts[to], pwgts[me], pwgts[me]);
         INC_DEC(pvol[to], pvol[me], pvol[me]);
         psurf[to] = psurf[to] + psurf[me] - 2.0*degrees[jbest];
         psurf[me] = 0;

         /* If we moved from/to the pmax subdomain find the new one! */
         if (me == pmax || to == pmax) 
           for (pmax=-1, i = 0; i < nparts; i++)
              if (psurf[i] != 0) {
                new = ARATIO(dim, psurf[i], pvol[i]);
                if (pmax == -1 || new > maxar) {
                  maxar = new;
                  pmax = i;
                }
              }

         IFSET(ctrl->dbglvl, DBG_MERGE,
               printf("Ndeg=%d Merge %d with %d jbest2=%d best=%f pmax=%d maxar=%f\n", 
                      degrees, me, to, jbest2, best, pmax, maxar));

       }
     }

     for (j = 0; j < ndegrees; j++)
        phtable[ptarget[j]] = -1;
  }

  /* Remove empty elements */
  for (i = 0; i < nvtxs; i++) {
     elpairs[i].key = where[i];
     elpairs[i].val = i;
  }
  idxkeysort(nvtxs, elpairs);

  for (last = 0, where[elpairs[0].val] = last, i = 1; i < nvtxs; i++) {
     if (elpairs[i].key > elpairs[i-1].key)
       last++;
     where[elpairs[i].val] = last;
  }
  ctrl->nparts  = last+1;

  graph->nmoves = nmoves;

  IMfree(&ptr, &ind, &psurf, &pvol, &pwgts, &elpairs, &ptarget, &degrees, &phtable, LTERM);
}


/*************************************************************************
* This function attracts from adjacent elements to fix small fused elements 
* The objective is to directly minimize the sum of the aspect ratios
**************************************************************************/
void Contribute_ARatio(CtrlType *ctrl, GraphType *graph)
{
  int i, j, l, me, to, el, nsels, nrels, ndegrees;
  int jbest;
  int dim, nvtxs, nparts, minsize, maxsize;
  idxtype *xadj, *vwgt, *adjncy, *where, *pwgts;
  idxtype *phtable, *ptarget;
  idxtype *ind, *ptr;
  realtype old, new, OldFromAR, OldToAR, NewFromAR, NewToAR;
  realtype best, id, ed;
  realtype *vvol, *vsurf, *adjwgt, *pvol, *psurf, *degrees;
  idxKeyValueType *elpairs;

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  vwgt   = graph->vwgt;
  vvol   = graph->vvol;
  vsurf  = graph->vsurf;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where  = graph->where;

  dim     = ctrl->dim;
  minsize = ctrl->minsize;
  maxsize = ctrl->maxsize;
  nparts  = ctrl->nparts;

  /* Create a part-to-vertex mapping; compute pwgts, psurf, pvol */
  ptr   = idxsmalloc(nparts+1, 0, "Contribute_ARatio: ptr");
  ind   = idxmalloc(nvtxs, "Contribute_ARatio: ind");
  pwgts = idxsmalloc(nparts, 0, "Contribute_ARatio: pwgts");
  pvol  = realsmalloc(nparts, 0.0, "Contribute_ARatio: pvol");
  psurf = realsmalloc(nparts, 0.0, "Contribute_ARatio: psurf");

  for (i = 0; i < nvtxs; i++) {
     me = where[i];
     ptr[me]++;
     pwgts[me] += vwgt[i];
     pvol[me] += vvol[i];
     psurf[me] += vsurf[i];
     for (j = xadj[i]; j < xadj[i+1]; j++)
        if (where[adjncy[j]] != me)
          psurf[me] += adjwgt[j];
  }

  MAKECSR(i, nparts, ptr);
  for (i = 0; i < nvtxs; i++) {
     me = where[i];
     ind[ptr[me]] = i;
     ptr[me]++;
  }

  for (i = nparts; i > 0; i--)
     ptr[i] = ptr[i-1];
  ptr[0] = 0;

  /* Find elements that are > minsize and sort them according to pwgts */
  /* nrels = number of regular elements */
  elpairs =(idxKeyValueType *) IMmalloc(nvtxs*sizeof(idxKeyValueType), "Contribute_ARatio:elpairs");
  for (nsels=nrels=i=0; i < nparts; i++) {
     if (pwgts[i] > minsize) {
       elpairs[nrels].key = pwgts[i];
       elpairs[nrels++].val = i;
     }
     else if (pwgts[i] < minsize) 
       nsels++;
  }
  idxkeysort(nrels, elpairs);

  IFSET(ctrl->dbglvl, DBG_CONTR,
        printf("===== nsels = %d ===== nrels = %d ===== nparts = %d =====\n",
               nsels, nrels, nparts));

  if (nsels == 0) {
    graph->nmoves=-1;
    IMfree(&ptr, &ind, &psurf, &pvol, &pwgts, &elpairs, LTERM);
    return;
  }

  graph->nmoves=0;

  ptarget = idxmalloc(nparts, "FusedElementGraph: ptarget");
  degrees = realsmalloc(nparts, 0, "FusedElementGraph: degrees");
  phtable = idxsmalloc(nparts, -1, "FusedElementGraph: phtable"); 

  /* Determine the connectivity of the regular fused elements that afford to lose */
  for (el = nrels-1; el >= 0; el--) {
     for (me=elpairs[el].val, l = ptr[me]; l < ptr[me+1]; l++) {
	i = ind[l];
        if (pwgts[me]-vwgt[i] < minsize)
          continue;

        for (id=ed=0.0, ndegrees=0, j = xadj[i]; j < xadj[i+1]; j++) {
           to = where[adjncy[j]];

	   if (to == me)
             id += adjwgt[j];
	   else
             ed += adjwgt[j];
       
           if (to == me || pwgts[to] >= minsize || pwgts[to]+vwgt[i] > maxsize)
             continue;

           if (phtable[to] == -1) {        /* connection is not created yet */
             ptarget[ndegrees] = to;
             degrees[ndegrees] = adjwgt[j]; 
             phtable[to] = ndegrees++;
           }
           else                            /* connection is already there */
             degrees[phtable[to]] += adjwgt[j];
        }

        /* Determine which of the ndegrees moves is the best */
        if (ndegrees > 0) {
          j = 0;
          to = ptarget[j];

	  OldFromAR = ARATIO(dim, psurf[me], pvol[me]);
	  OldToAR   = ARATIO(dim, psurf[to], pvol[to]);
          NewFromAR = ARATIO(dim, psurf[me]+id-ed-vsurf[i], pvol[me]-vvol[i]);
          NewToAR   = ARATIO(dim, psurf[to]+ed+id-2.0*degrees[j]+vsurf[i], pvol[to]+vvol[i]);

	  old = OldFromAR + OldToAR;
          new = NewFromAR + NewToAR;

	  jbest = j;
	  best = old - new;

          for (j = 1; j < ndegrees; j++) {
             to = ptarget[j];

	     OldToAR = ARATIO(dim, psurf[to], pvol[to]);
             NewToAR = ARATIO(dim, psurf[to]+ed+id-2.0*degrees[j]+vsurf[i], pvol[to]+vvol[i]);
	    
	     old = OldFromAR + OldToAR;
             new = NewFromAR + NewToAR;
   
	     /* Check if it increases the aspect ratio */
             if (best < old - new) {
               jbest = j;
               best = old - new;
             }
          }
   
          /* Contribute now */
          to = ptarget[jbest];
	  where[i] = to;

	  IFSET(ctrl->dbglvl, DBG_CONTR,
                printf("Ndeg=%d Move %d from %d to %d best=%f\n", ndegrees, i,
                       me, to, best));
   
          INC_DEC(pwgts[to], pwgts[me], vwgt[i]);
          INC_DEC(pvol[to], pvol[me], vvol[i]);
          psurf[me] = psurf[me] + id - ed - vsurf[i];
          psurf[to] = psurf[to] + id + ed - 2.0*degrees[jbest] + vsurf[i];
        }

        for (j = 0; j < ndegrees; j++)
           phtable[ptarget[j]] = -1;
     }
  }

  IMfree(&ptr, &ind, &psurf, &pvol, &pwgts, &elpairs, &ptarget, &degrees, &phtable, LTERM);
}


/****************************************************************************
* This function attracts from adjacent elements to fix small fused elements 
* The objective is to directly minimize the weighted sum of the aspect ratios
*****************************************************************************/
void Contribute_WeightARatio(CtrlType *ctrl, GraphType *graph)
{
  int i, j, l, me, to, el, nsels, nrels, ndegrees;
  int jbest;
  int dim, nvtxs, nparts, minsize, maxsize;
  idxtype *xadj, *vwgt, *adjncy, *where, *pwgts;
  idxtype *phtable, *ptarget;
  idxtype *ind, *ptr;
  realtype old, new, OldFromAR, OldToAR, NewFromAR, NewToAR;
  realtype best, id, ed;
  realtype *vvol, *vsurf, *adjwgt, *pvol, *psurf, *degrees;
  idxKeyValueType *elpairs;

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  vwgt   = graph->vwgt;
  vvol   = graph->vvol;
  vsurf  = graph->vsurf;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where  = graph->where;

  dim     = ctrl->dim;
  minsize = ctrl->minsize;
  maxsize = ctrl->maxsize;
  nparts  = ctrl->nparts;

  /* Create a part-to-vertex mapping; compute pwgts, psurf, pvol */
  ptr   = idxsmalloc(nparts+1, 0, "Contribute_ARatio: ptr");
  ind   = idxmalloc(nvtxs, "Contribute_ARatio: ind");
  pwgts = idxsmalloc(nparts, 0, "Contribute_ARatio: pwgts");
  pvol  = realsmalloc(nparts, 0.0, "Contribute_ARatio: pvol");
  psurf = realsmalloc(nparts, 0.0, "Contribute_ARatio: psurf");

  for (i = 0; i < nvtxs; i++) {
     me = where[i];
     ptr[me]++;
     pwgts[me] += vwgt[i];
     pvol[me] += vvol[i];
     psurf[me] += vsurf[i];
     for (j = xadj[i]; j < xadj[i+1]; j++)
        if (where[adjncy[j]] != me)
          psurf[me] += adjwgt[j];
  }

  MAKECSR(i, nparts, ptr);
  for (i = 0; i < nvtxs; i++) {
     me = where[i];
     ind[ptr[me]] = i;
     ptr[me]++;
  }

  for (i = nparts; i > 0; i--)
     ptr[i] = ptr[i-1];
  ptr[0] = 0;

  /* Find elements that are > minsize and sort them according to pwgts */
  /* nrels = number of regular elements */
  elpairs =(idxKeyValueType *) IMmalloc(nparts*sizeof(idxKeyValueType), "Contribute_ARatio:elpairs");
  for (nsels=nrels=i=0; i < nparts; i++) {
     if (pwgts[i] > minsize) {
       elpairs[nrels].key = pwgts[i];
       elpairs[nrels++].val = i;
     }
     else if (pwgts[i] < minsize) 
       nsels++;
  }
  idxkeysort(nrels, elpairs);

  IFSET(ctrl->dbglvl, DBG_CONTR,
        printf("===== nsels = %d ===== nrels = %d ===== nparts = %d =====\n",
               nsels, nrels, nparts));

  if (nsels == 0) {
    graph->nmoves=-1;
    IMfree(&ptr, &ind, &psurf, &pvol, &pwgts, &elpairs, LTERM);
    return;
  }

  graph->nmoves=0;

  ptarget = idxmalloc(nparts, "FusedElementGraph: ptarget");
  degrees = realsmalloc(nparts, 0, "FusedElementGraph: degrees");
  phtable = idxsmalloc(nparts, -1, "FusedElementGraph: phtable"); 

  /* Determine the connectivity of the regular fused elements that afford to lose */
  for (el = nrels-1; el >= 0; el--) {
     for (me=elpairs[el].val, l = ptr[me]; l < ptr[me+1]; l++) {
	i = ind[l];
        if (pwgts[me]-vwgt[i] < minsize)
          continue;

        for (id=ed=0.0, ndegrees=0, j = xadj[i]; j < xadj[i+1]; j++) {
           to = where[adjncy[j]];

	   if (to == me)
             id += adjwgt[j];
	   else
             ed += adjwgt[j];
       
           if (to == me || pwgts[to] >= minsize || pwgts[to]+vwgt[i] > maxsize)
             continue;

           if (phtable[to] == -1) {        /* connection is not created yet */
             ptarget[ndegrees] = to;
             degrees[ndegrees] = adjwgt[j]; 
             phtable[to] = ndegrees++;
           }
           else                            /* connection is already there */
             degrees[phtable[to]] += adjwgt[j];
        }

        /* Determine which of the ndegrees moves is the best */
        if (ndegrees > 0) {
          j = 0;
          to = ptarget[j];

	  OldFromAR = ARATIO(dim, psurf[me], pvol[me]) * pwgts[me];
	  OldToAR   = ARATIO(dim, psurf[to], pvol[to]);
          NewFromAR = ARATIO(dim, psurf[me]+id-ed-vsurf[i], pvol[me]-vvol[i]);
          NewToAR   = ARATIO(dim, psurf[to]+ed+id-2.0*degrees[j]+vsurf[i], pvol[to]+vvol[i]);

	  old = OldFromAR + OldToAR*pwgts[to];
          new = NewFromAR*(pwgts[me]-vwgt[i]) + NewToAR*(pwgts[to]+vwgt[i]);

	  jbest = j;
	  best = old - new;

          for (j = 1; j < ndegrees; j++) {
             to = ptarget[j];

	     OldToAR = ARATIO(dim, psurf[to], pvol[to]);
             NewToAR = ARATIO(dim, psurf[to]+ed+id-2.0*degrees[j]+vsurf[i], pvol[to]+vvol[i]);
	    
	     old = OldFromAR + OldToAR*pwgts[to];
             new = NewFromAR*(pwgts[me]-vwgt[i]) + NewToAR*(pwgts[to]+vwgt[i]);
   
	     /* Check if it increases the aspect ratio */
             if (best < old - new) {
               jbest = j;
               best = old - new;
             }
          }
   
          /* Contribute now */
          to = ptarget[jbest];
	  where[i] = to;

	  IFSET(ctrl->dbglvl, DBG_CONTR,
                printf("Ndeg=%d Move %d from %d to %d best=%f\n", ndegrees, i,
                me, to, best));
   
          INC_DEC(pwgts[to], pwgts[me], vwgt[i]);
          INC_DEC(pvol[to], pvol[me], vvol[i]);
          psurf[me] = psurf[me] + id - ed - vsurf[i];
          psurf[to] = psurf[to] + id + ed - 2.0*degrees[jbest] + vsurf[i];
   
        }

        for (j = 0; j < ndegrees; j++)
           phtable[ptarget[j]] = -1;
     }
  }

  IMfree(&ptr, &ind, &psurf, &pvol, &pwgts, &elpairs, &ptarget, &degrees, &phtable, LTERM);
}


/*************************************************************************
* This function attracts from adjacent elements to fix small fused elements 
* The objective is to directly minimize the maximum aspect ratio.
**************************************************************************/
void Contribute_MinMaxARatio(CtrlType *ctrl, GraphType *graph)
{
  int i, j, l, me, to, el, nsels, nrels, ndegrees, pmax;
  int jbest;
  int dim, nvtxs, nparts, minsize, maxsize;
  idxtype *xadj, *vwgt, *adjncy, *where, *pwgts;
  idxtype *phtable, *ptarget;
  idxtype *ind, *ptr;
  realtype new, maxar, NewFromAR, NewToAR;
  realtype best, id, ed;
  realtype *vvol, *vsurf, *adjwgt, *pvol, *psurf, *degrees;
  idxKeyValueType *elpairs;

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  vwgt   = graph->vwgt;
  vvol   = graph->vvol;
  vsurf  = graph->vsurf;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where  = graph->where;

  dim     = ctrl->dim;
  minsize = ctrl->minsize;
  maxsize = ctrl->maxsize;
  nparts  = ctrl->nparts;

  /* Create a part-to-vertex mapping; compute pwgts, psurf, pvol */
  ptr   = idxsmalloc(nparts+1, 0, "Contribute_MinMaxARatio: ptr");
  ind   = idxmalloc(nvtxs, "Contribute_MinMaxARatio: ind");
  pwgts = idxsmalloc(nparts, 0, "Contribute_MinMaxARatio: pwgts");
  pvol  = realsmalloc(nparts, 0.0, "Contribute_MinMaxARatio: pvol");
  psurf = realsmalloc(nparts, 0.0, "Contribute_MinMaxARatio: psurf");

  for (i = 0; i < nvtxs; i++) {
     me = where[i];
     ptr[me]++;
     pwgts[me] += vwgt[i];
     pvol[me] += vvol[i];
     psurf[me] += vsurf[i];
     for (j = xadj[i]; j < xadj[i+1]; j++)
        if (where[adjncy[j]] != me)
          psurf[me] += adjwgt[j];
  }

  MAKECSR(i, nparts, ptr);
  for (i = 0; i < nvtxs; i++) {
     me = where[i];
     ind[ptr[me]] = i;
     ptr[me]++;
  }

  for (i = nparts; i > 0; i--)
     ptr[i] = ptr[i-1];
  ptr[0] = 0;

  /* Find elements that are > minsize and sort them according to pwgts */
  /* nrels = number of regular elements */
  elpairs =(idxKeyValueType *) IMmalloc(nparts*sizeof(idxKeyValueType), "Merge:elpairs");
  for (nsels=nrels=i=0; i < nparts; i++) {
     if (pwgts[i] > minsize) {
       elpairs[nrels].key = pwgts[i];
       elpairs[nrels++].val = i;
     }
     else if (pwgts[i] < minsize) 
       nsels++;
  }
  idxkeysort(nrels, elpairs);

  IFSET(ctrl->dbglvl, DBG_CONTR,
        printf("===== nsels = %d ===== nrels = %d ===== nparts = %d =====\n",
               nsels, nrels, nparts));

  if (nsels == 0) {
    graph->nmoves=-1;
    IMfree(&ptr, &ind, &psurf, &pvol, &pwgts, &elpairs, LTERM);
    return;
  }

  graph->nmoves=0;

  ptarget = idxmalloc(nparts, "FusedElementGraph: ptarget");
  degrees = realsmalloc(nparts, 0, "FusedElementGraph: degrees");
  phtable = idxsmalloc(nparts, -1, "FusedElementGraph: phtable"); 

  /* Determine the domain that has the maximum aspect ratio */
  pmax = 0;
  maxar = ARATIO(dim, psurf[0], pvol[0]);
  for (i = 1; i < nparts; i++) {
     new = ARATIO(dim, psurf[i], pvol[i]);
     if (new > maxar) {
       maxar = new;
       pmax = i;
     }
  }
  IFSET(ctrl->dbglvl, DBG_CONTR,
        printf("===== pmax = %d ===== maxar = %f =====n",pmax, maxar));

  /* Determine the connectivity of the regular fused elements that afford to lose */
  for (el = nrels-1; el >= 0; el--) {
     for (me=elpairs[el].val, l = ptr[me]; l < ptr[me+1]; l++) {
	i = ind[l];
        if (pwgts[me]-vwgt[i] < minsize)
          continue;

        for (id=ed=0.0, ndegrees=0, j = xadj[i]; j < xadj[i+1]; j++) {
           to = where[adjncy[j]];

	   if (to == me)
             id += adjwgt[j];
	   else
             ed += adjwgt[j];
       
           if (to == me || pwgts[to] >= minsize || pwgts[to]+vwgt[i] > maxsize)
             continue;

           if (phtable[to] == -1) {        /* connection is not created yet */
             ptarget[ndegrees] = to;
             degrees[ndegrees] = adjwgt[j]; 
             phtable[to] = ndegrees++;
           }
           else                            /* connection is already there */
             degrees[phtable[to]] += adjwgt[j];
        }

        /* Determine which of the ndegrees moves is the best */
        if (ndegrees > 0) {
          j = 0;
          to = ptarget[j];

          NewFromAR = ARATIO(dim, psurf[me]+id-ed-vsurf[i], pvol[me]-vvol[i]);
          NewToAR   = ARATIO(dim, psurf[to]+ed+id-2.0*degrees[j]+vsurf[i], pvol[to]+vvol[i]);

	  jbest = j;
	  best = amax(NewFromAR, NewToAR);

          for (j = 1; j < ndegrees; j++) {
             to = ptarget[j];

             NewToAR = ARATIO(dim, psurf[to]+ed+id-2.0*degrees[j]+vsurf[i], pvol[to]+vvol[i]);
	     new = amax(NewFromAR, NewToAR);
   
	     /* Check if it increases the aspect ratio */
             if (new < best) {
               jbest = j;
               best = new;
             }
          }
   
          /* Merge now */
          to = ptarget[jbest];
	  where[i] = to;

          INC_DEC(pwgts[to], pwgts[me], vwgt[i]);
          INC_DEC(pvol[to], pvol[me], vvol[i]);
          psurf[me] = psurf[me] + id - ed - vsurf[i];
          psurf[to] = psurf[to] + id + ed - 2.0*degrees[jbest] + vsurf[i];
   
          /* find the new maximum aspect ratio */
          if ( best > maxar || me == pmax || to == pmax) {
            pmax = 0;
            maxar = ARATIO(dim, psurf[0], pvol[0]);
            for (i = 1; i < nparts; i++)
               if ((new = ARATIO(dim, psurf[i], pvol[i])) > maxar) {
                 maxar = new;
                 pmax = i;
               }
          }

	  IFSET(ctrl->dbglvl, DBG_CONTR, 
                printf("Ndeg=%d Move %d from %d to %d jbest=%d best=%f pmax=%d maxar=%f\n",
                       ndegrees, i, me, to, jbest, best, pmax, maxar));
   
        }

        for (j = 0; j < ndegrees; j++)
           phtable[ptarget[j]] = -1;
     }
  }

  IMfree(&ptr, &ind, &psurf, &pvol, &pwgts, &elpairs, &ptarget, &degrees, &phtable, LTERM);
}


/*************************************************************************
* This function attracts from adjacent elements to fix small fused elements 
* The objective is to directly minimize the multiple objective aspect ratio
* In this case RType = REFINE_MINMAXAR and RType =  REFINE_WAR
**************************************************************************/
void Contribute_MultiObj(CtrlType *ctrl, GraphType *graph)
{
  int i, j, l, me, to, el, nsels, nrels, ndegrees, pmax;
  int jbest1, jbest2, jbest;
  int dim, nvtxs, nparts, minsize, maxsize;
  idxtype *xadj, *vwgt, *adjncy, *where, *pwgts;
  idxtype *phtable, *ptarget;
  idxtype *ind, *ptr;
  realtype old, new, maxar, OldFromAR, OldToAR, NewFromAR, NewToAR;
  realtype best1, best2, best, id, ed;
  realtype *vvol, *vsurf, *adjwgt, *pvol, *psurf, *degrees;
  idxKeyValueType *elpairs;

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  vwgt   = graph->vwgt;
  vvol   = graph->vvol;
  vsurf  = graph->vsurf;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where  = graph->where;

  dim     = ctrl->dim;
  minsize = ctrl->minsize;
  maxsize = ctrl->maxsize;
  nparts  = ctrl->nparts;

  /* Create a part-to-vertex mapping; compute pwgts, psurf, pvol */
  ptr   = idxsmalloc(nparts+1, 0, "Merge: ptr");
  ind   = idxmalloc(nvtxs, "Merge: ind");
  pwgts = idxsmalloc(nparts, 0, "Merge: pwgts");
  pvol  = realsmalloc(nparts, 0.0, "Merge: pvol");
  psurf = realsmalloc(nparts, 0.0, "Merge: psurf");

  for (i = 0; i < nvtxs; i++) {
     me = where[i];
     ptr[me]++;
     pwgts[me] += vwgt[i];
     pvol[me] += vvol[i];
     psurf[me] += vsurf[i];
     for (j = xadj[i]; j < xadj[i+1]; j++)
        if (where[adjncy[j]] != me)
          psurf[me] += adjwgt[j];
  }

  MAKECSR(i, nparts, ptr);
  for (i = 0; i < nvtxs; i++) {
     me = where[i];
     ind[ptr[me]] = i;
     ptr[me]++;
  }

  for (i = nparts; i > 0; i--)
     ptr[i] = ptr[i-1];
  ptr[0] = 0;

  /* Find elements that are > minsize and sort them according to pwgts */
  /* nrels = number of regular elements */
  elpairs =(idxKeyValueType *) IMmalloc(nparts*sizeof(idxKeyValueType), "Merge:elpairs");
  for (nsels=nrels=i=0; i < nparts; i++) {
     if (pwgts[i] > minsize) {
       elpairs[nrels].key = pwgts[i];
       elpairs[nrels++].val = i;
     }
     else if (pwgts[i] < minsize) 
       nsels++;
  }
  idxkeysort(nrels, elpairs);

  IFSET(ctrl->dbglvl, DBG_CONTR,
        printf("===== nsels = %d ===== nrels = %d ===== nparts = %d =====\n",
               nsels, nrels, nparts));

  if (nsels == 0) {
    graph->nmoves=-1;
    IMfree(&ptr, &ind, &psurf, &pvol, &pwgts, &elpairs, LTERM);
    return;
  }

  graph->nmoves=0;

  ptarget = idxmalloc(nparts, "FusedElementGraph: ptarget");
  degrees = realsmalloc(nparts, 0, "FusedElementGraph: degrees");
  phtable = idxsmalloc(nparts, -1, "FusedElementGraph: phtable"); 

  /* Determine the domain that has the maximum aspect ratio */
  pmax = 0;
  maxar = ARATIO(dim, psurf[0], pvol[0]);
  for (i = 1; i < nparts; i++) {
     new = ARATIO(dim, psurf[i], pvol[i]);
     if (new > maxar) {
       maxar = new;
       pmax = i;
     }
  }

  /* Determine the connectivity of the regular fused elements that afford to lose */
  for (el = nrels-1; el >= 0; el--) {
     for (me=elpairs[el].val, l = ptr[me]; l < ptr[me+1]; l++) {
	i = ind[l];
        if (pwgts[me]-vwgt[i] < minsize)
          continue;

        for (id=ed=0.0, ndegrees=0, j = xadj[i]; j < xadj[i+1]; j++) {
           to = where[adjncy[j]];

	   if (to == me)
             id += adjwgt[j];
	   else
             ed += adjwgt[j];
       
           if (to == me || pwgts[to] >= minsize || pwgts[to]+vwgt[i] > maxsize)
/* if (to == me || pwgts[to] >= pwgts[me] || pwgts[to]+vwgt[i] > maxsize) */
             continue;

           if (phtable[to] == -1) {        /* connection is not created yet */
             ptarget[ndegrees] = to;
             degrees[ndegrees] = adjwgt[j]; 
             phtable[to] = ndegrees++;
           }
           else                            /* connection is already there */
             degrees[phtable[to]] += adjwgt[j];
        }

        /* Determine which of the ndegrees moves is the best */
        if (ndegrees > 0) {
          j = 0;
          to = ptarget[j];

          OldFromAR = ARATIO(dim, psurf[me], pvol[me]) * pwgts[me];
          NewFromAR = ARATIO(dim, psurf[me]+id-ed-vsurf[i], pvol[me]-vvol[i]);
          NewToAR   = ARATIO(dim, psurf[to]+ed+id-2.0*degrees[j]+vsurf[i], pvol[to]+vvol[i]);

	  jbest1 = j;
	  best1 = amax(NewFromAR, NewToAR);

	  if (best1 <= maxar) {
            OldToAR   = ARATIO(dim, psurf[to], pvol[to]);

	    old = OldFromAR + OldToAR*pwgts[to];
	    new = NewFromAR*(pwgts[me]-vwgt[i]) + NewToAR*(pwgts[to]+vwgt[i]);

            jbest2 = jbest1;
            best2 = old - new;
          }
          else
            jbest2 = -1;

          for (j = 1; j < ndegrees; j++) {
             to = ptarget[j];

             NewToAR = ARATIO(dim, psurf[to]+ed+id-2.0*degrees[j]+vsurf[i], pvol[to]+vvol[i]);
	     new = amax(NewFromAR, NewToAR);
   
	     /* Check if it increases the aspect ratio */
             if (new < best1) {
               jbest1 = j;
               best1 = new;
             }
   
             /* If first objective is OK, check second one */
             if (new <= maxar) {
               OldToAR   = ARATIO(dim, psurf[to], pvol[to]);
   
	       old = OldFromAR + OldToAR*pwgts[to];
               new = NewFromAR*(pwgts[me]-vwgt[i]) + NewToAR*(pwgts[to]+vwgt[i]);

               if ( jbest2 == -1 || best2 < old-new ) {
                 jbest2 = j;
                 best2  = old-new;
               }
             }
      
          }
      
          if (jbest2 != -1) {
            jbest = jbest2;
            best = best2;
          }
          else {
            jbest = jbest1;
            best = best1;
          }
   
          /* Merge now */
          to = ptarget[jbest];
	  where[i] = to;

          INC_DEC(pwgts[to], pwgts[me], vwgt[i]);
          INC_DEC(pvol[to], pvol[me], vvol[i]);
          psurf[me] = psurf[me] + id - ed - vsurf[i];
          psurf[to] = psurf[to] + id + ed - 2.0*degrees[jbest] + vsurf[i];
   
          /* find the new maximum aspect ratio */
          if (jbest2 == -1 || me == pmax || to == pmax) {
            pmax = 0;
            maxar = ARATIO(dim, psurf[0], pvol[0]);
            for (i = 1; i < nparts; i++)
               if ((new = ARATIO(dim, psurf[i], pvol[i])) > maxar) {
                 maxar = new;
                 pmax = i;
               }
          }

          IFSET(ctrl->dbglvl, DBG_CONTR, 
                printf("Ndeg=%d Move %d from %d to %d jbest2=%d best=%f pmax=%d maxar=%f\n",
                       ndegrees, i, me, to, jbest2, best, pmax, maxar));
        }

        for (j = 0; j < ndegrees; j++)
           phtable[ptarget[j]] = -1;
     }
  }

  IMfree(&ptr, &ind, &psurf, &pvol, &pwgts, &elpairs, &ptarget, &degrees, &phtable, LTERM);
}
