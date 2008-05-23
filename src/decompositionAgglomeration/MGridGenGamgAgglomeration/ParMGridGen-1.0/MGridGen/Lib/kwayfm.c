/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * kwayfm.c
 *
 * This file contains code that implements the multilevel k-way refinement
 *
 * George Irene
 */

#include "mgridgen.h"


/*************************************************************************
* This function performs k-way refinement, whose objective is to directly
* minimize the sum of the aspect ratios
**************************************************************************/
void Random_KWayARatioRefine(CtrlType *ctrl, GraphType *graph, int npasses)
{
  int i, ii, j, dim, nparts, pass, nvtxs, nmoves, ndegrees, pmax; 
  int from, to, jbest;
  idxtype *xadj, *vwgt, *adjncy, *where, *pwgts, *perm, *phtable, *ptarget;
  realtype old, new, best, id, ed, maxar;
  realtype *vvol, *vsurf, *adjwgt, *pvol, *psurf, *degrees;

  nvtxs     = graph->nvtxs;
  xadj      = graph->xadj;
  vwgt      = graph->vwgt;
  vvol      = graph->vvol;
  vsurf     = graph->vsurf;
  adjncy    = graph->adjncy;
  adjwgt    = graph->adjwgt;

  where  = graph->where;
  pwgts  = graph->pwgts;
  pvol   = graph->pvol;
  psurf  = graph->psurf;
  
  dim	   = ctrl->dim;
  nparts   = ctrl->nparts;
  degrees  = realmalloc(nparts, "degrees");
  phtable  = idxsmalloc(nparts, -1, "phtable");
  ptarget  = idxsmalloc(nparts, -1, "ptarget");
  perm     = idxmalloc(nvtxs, "perm");

  /* Determine the domain that has the maximum aspect ratio */
  pmax = 0;
  maxar = ARATIO(dim, psurf[0], pvol[0]);
  for (i=1; i<nparts; i++) {
     new = ARATIO(dim, psurf[i], pvol[i]);
     if (new > maxar) {
       maxar = new;
       pmax = i;
     }
  }

  IFSET(ctrl->dbglvl, DBG_REFINE,
     printf("Partitions: [%3d %3d]-[%3d %3d]. MaxRatio: [%4d, %e], Ratio: %e\n",
             pwgts[iamin(nparts, pwgts)], pwgts[iamax(nparts, pwgts)], 
             ctrl->minsize, ctrl->maxsize, pmax, maxar, graph->minratio));

  RandomPermute(nvtxs, perm, 1);
  for (pass=0; pass<npasses; pass++) {
     RandomPermute(nvtxs, perm, 0);
     RandomPermute(nvtxs, perm, 0);
     for (nmoves=ii=0; ii<nvtxs; ii++) {
        i = perm[ii];
        from = where[i];

        if (pwgts[from] - vwgt[i] < ctrl->minsize)
          continue;

        /* Determine the connectivity of the 'i' vertex */
        for (id=ed=0.0, ndegrees=0, j=xadj[i]; j<xadj[i+1]; j++) {
           to = where[adjncy[j]];

           if (to == from)
             id += adjwgt[j];
           else
             ed += adjwgt[j];

           if (to != from && pwgts[to]+vwgt[i] <= ctrl->maxsize) {
             if (phtable[to] == -1) {
               degrees[ndegrees] = adjwgt[j];
               ptarget[ndegrees] = to;
               phtable[to] = ndegrees++;
             }
             else 
               degrees[phtable[to]] += adjwgt[j];
           }
        }

        /* Determine which of the ndegrees moves is the best */
        for (best=0.1, jbest=-1, j=0; j<ndegrees; j++) {
           to = ptarget[j];
           old = ARATIO(dim, psurf[from], pvol[from]) + ARATIO(dim, psurf[to], pvol[to]);
           new = ARATIO(dim, psurf[from]+id-ed-vsurf[i], pvol[from]-vvol[i]) +
                 ARATIO(dim, psurf[to]+ed+id-2.0*degrees[j]+vsurf[i], pvol[to]+vvol[i]);

           if (best < old-new) {
             best = old-new;
             jbest = j;
           }
        }

        if (jbest != -1) {
          to = ptarget[jbest];
          where[i] = to;
          INC_DEC(pwgts[to], pwgts[from], vwgt[i]); 
          INC_DEC(pvol[to], pvol[from], vvol[i]);
          psurf[from] = psurf[from] + id - ed - vsurf[i];
          psurf[to]   = psurf[to] + id + ed - 2.0*degrees[jbest] + vsurf[i];
          graph->minratio = ComputeFunction(ctrl->RType, ctrl, graph);
          nmoves++;

          IFSET(ctrl->dbglvl, DBG_MOVEINFO,
            printf("Moving %6d from %3d to %3d. Gain: %8.6f. MinRatio: %e [%e]\n"
                   , i, from, to, best, graph->minratio, vsurf[i]));

          /* CheckParams(ctrl, graph); */
        }

        for (j=0; j<ndegrees; j++) 
           phtable[ptarget[j]] = -1;
     }

     IFSET(ctrl->dbglvl, DBG_REFINE,
           printf("\t[%6d %6d], Nmoves: %5d, MinRatio: %e\n",
           pwgts[iamin(nparts, pwgts)], pwgts[iamax(nparts, pwgts)], nmoves,
           graph->minratio));

     if (nmoves == 0)
       break;
  }

  pmax = 0;
  maxar = ARATIO(dim, psurf[0], pvol[0]);
  for (i=1; i<nparts; i++) {
     new = ARATIO(dim, psurf[i], pvol[i]);
     if (new > maxar) {
       maxar = new;
       pmax = i;
     }
  }

  graph->nmoves = nmoves;
  IFSET(ctrl->dbglvl, DBG_REFINE, printf("FinalMax: %d %e\n", pmax, maxar));

  IMfree(&perm, &phtable, &degrees, &ptarget, LTERM);
}


/*************************************************************************
* This function performs k-way refinement, whose objective is to directly
* minimize the weighted sum of the aspect ratios
**************************************************************************/
void Random_KWayWeightARatioRefine(CtrlType *ctrl, GraphType *graph, int npasses)
{
  int i, ii, j, dim, nparts, pass, nvtxs, nmoves, ndegrees, pmax; 
  int from, to, jbest;
  idxtype *xadj, *vwgt, *adjncy, *where, *pwgts, *perm, *phtable, *ptarget;
  realtype old, new, best, id, ed, maxar;
  realtype *vvol, *vsurf, *adjwgt, *pvol, *psurf, *degrees;

  nvtxs     = graph->nvtxs;
  xadj      = graph->xadj;
  vwgt      = graph->vwgt;
  vvol      = graph->vvol;
  vsurf     = graph->vsurf;
  adjncy    = graph->adjncy;
  adjwgt    = graph->adjwgt;

  where  = graph->where;
  pwgts  = graph->pwgts;
  pvol   = graph->pvol;
  psurf  = graph->psurf;
  
  dim      = ctrl->dim;
  nparts   = ctrl->nparts;
  degrees  = realmalloc(nparts, "degrees");
  phtable  = idxsmalloc(nparts, -1, "phtable");
  ptarget  = idxsmalloc(nparts, -1, "ptarget");
  perm     = idxmalloc(nvtxs, "perm");

  /* Determine the domain that has the maximum aspect ratio */
  pmax = 0;
  maxar = ARATIO(dim, psurf[0], pvol[0]);
  for (i=1; i<nparts; i++) {
     new = ARATIO(dim, psurf[i], pvol[i]);
     if (new > maxar) {
       maxar = new;
       pmax = i;
     }
  }

  IFSET(ctrl->dbglvl, DBG_REFINE,
        printf("Partitions: [%3d %3d]-[%3d %3d].  MaxRatio: [%4d, %e], Ratio: %e\n",
         pwgts[iamin(nparts, pwgts)], pwgts[iamax(nparts, pwgts)],
         ctrl->minsize, ctrl->maxsize, pmax, maxar, graph->minratio));

  RandomPermute(nvtxs, perm, 1);
  for (pass=0; pass<npasses; pass++) {
     RandomPermute(nvtxs, perm, 0);
     RandomPermute(nvtxs, perm, 0);
     for (nmoves=ii=0; ii<nvtxs; ii++) {
        i = perm[ii];
        from = where[i];

        if (pwgts[from] - vwgt[i] < ctrl->minsize)
          continue;

        /* Determine the connectivity of the 'i' vertex */
        for (id=ed=0.0, ndegrees=0, j=xadj[i]; j<xadj[i+1]; j++) {
           to = where[adjncy[j]];

           if (to == from)
             id += adjwgt[j];
           else
             ed += adjwgt[j];
 
           if (to != from && pwgts[to]+vwgt[i] <= ctrl->maxsize) {
             if (phtable[to] == -1) {
               degrees[ndegrees] = adjwgt[j];
               ptarget[ndegrees] = to;
               phtable[to] = ndegrees++;
             }
             else 
               degrees[phtable[to]] += adjwgt[j];
           }
        }

        /* Determine which of the ndegrees moves is the best */
        for (best=.1, jbest=-1, j=0; j<ndegrees; j++) {
           to = ptarget[j];
           old = pwgts[from]*ARATIO(dim, psurf[from], pvol[from]) + 
                 pwgts[to]*ARATIO(dim, psurf[to], pvol[to]);
           new = (pwgts[from]-vwgt[i])*ARATIO(dim, psurf[from]+id-ed-vsurf[i], pvol[from]-vvol[i]) +
                (pwgts[to]+vwgt[i])*ARATIO(dim, psurf[to]+ed+id-2.0*degrees[j]+vsurf[i], pvol[to]+vvol[i]);
           if (best < old-new) {
             best = old-new;
             jbest = j;
           }
        }

        if (jbest != -1) {
          to = ptarget[jbest];
          where[i] = to;
          INC_DEC(pwgts[to], pwgts[from], vwgt[i]); 
          INC_DEC(pvol[to], pvol[from], vvol[i]);
          psurf[from] = psurf[from] + id - ed - vsurf[i];
          psurf[to]   = psurf[to] + id + ed - 2.0*degrees[jbest] + vsurf[i];
          graph->minratio = ComputeFunction(ctrl->RType, ctrl, graph);
          nmoves++;

          IFSET(ctrl->dbglvl, DBG_MOVEINFO,
            printf("\tMoving %6d from %3d to %3d. Gain: %4.2f. MinRatio: %e [%e]\n",
                   i, from, to, best, graph->minratio, vsurf[i]));
          /* CheckParams(ctrl, graph); */
        }

        for (j=0; j<ndegrees; j++) 
           phtable[ptarget[j]] = -1;
     }

     IFSET(ctrl->dbglvl, DBG_REFINE,
           printf("\t[%6d %6d], Nmoves: %5d, MinRatio: %e\n", pwgts[iamin(nparts, pwgts)],
                        pwgts[iamax(nparts, pwgts)], nmoves, graph->minratio));

     if (nmoves == 0)
       break;
  }

  pmax = 0;
  maxar = ARATIO(dim, psurf[0], pvol[0]);
  for (i=1; i<nparts; i++) {
     new = ARATIO(dim, psurf[i], pvol[i]);
     if (new > maxar) {
       maxar = new;
       pmax = i;
     } 
  }

  graph->nmoves = nmoves;
  IFSET(ctrl->dbglvl, DBG_REFINE, printf("FinalMax: %d %e\n", pmax, maxar));

  IMfree(&perm, &phtable, &degrees, &ptarget, LTERM);
}


/*************************************************************************
* This function performs k-way refinement, whose objective is to directly
* minimize the surface cut
**************************************************************************/
void Random_KWaySCutRefine(CtrlType *ctrl, GraphType *graph, int npasses)
{
  int i, ii, j, dim, nparts, pass, nvtxs, nmoves, ndegrees, pmax; 
  int from, to, jbest;
  idxtype *xadj, *vwgt, *adjncy, *where, *pwgts, *perm, *phtable, *ptarget;
  realtype old, new, best, id, ed, maxar;
  realtype *vvol, *vsurf, *adjwgt, *pvol, *psurf, *degrees;

  nvtxs     = graph->nvtxs;
  xadj      = graph->xadj;
  vwgt      = graph->vwgt;
  vvol      = graph->vvol;
  vsurf     = graph->vsurf;
  adjncy    = graph->adjncy;
  adjwgt    = graph->adjwgt;

  where  = graph->where;
  pwgts  = graph->pwgts;
  pvol   = graph->pvol;
  psurf  = graph->psurf;
  
  dim      = ctrl->dim;
  nparts   = ctrl->nparts;
  degrees  = realmalloc(nparts, "degrees");
  phtable  = idxsmalloc(nparts, -1, "phtable");
  ptarget  = idxsmalloc(nparts, -1, "ptarget");
  perm     = idxmalloc(nvtxs, "perm");


  /* Determine the domain that has the maximum aspect ratio */
  pmax = 0;
  maxar = ARATIO(dim, psurf[0], pvol[0]);
  for (i=1; i<nparts; i++) {
     new = ARATIO(dim, psurf[i], pvol[i]);
     if (new > maxar) {
       maxar = new;
       pmax = i;
     }
  }

  IFSET(ctrl->dbglvl, DBG_REFINE,
     printf("Partitions: [%3d %3d]-[%3d %3d]. MaxRatio: [%4d, %e], Ratio: %e\n",
             pwgts[iamin(nparts, pwgts)], pwgts[iamax(nparts, pwgts)], 
             ctrl->minsize, ctrl->maxsize, pmax, maxar, graph->minratio));

  RandomPermute(nvtxs, perm, 1);
  for (pass=0; pass<npasses; pass++) {
     RandomPermute(nvtxs, perm, 0);
     RandomPermute(nvtxs, perm, 0);
     for (nmoves=ii=0; ii<nvtxs; ii++) {
        i = perm[ii];
        from = where[i];

        if (pwgts[from] - vwgt[i] < ctrl->minsize)
          continue;

        /* Determine the connectivity of the 'i' vertex */
        for (id=ed=0.0, ndegrees=0, j=xadj[i]; j<xadj[i+1]; j++) {
           to = where[adjncy[j]];

           if (to == from)
             id += adjwgt[j];
           else
             ed += adjwgt[j];

           if (to != from && pwgts[to]+vwgt[i] <= ctrl->maxsize) {
             if (phtable[to] == -1) {
               degrees[ndegrees] = adjwgt[j];
               ptarget[ndegrees] = to;
               phtable[to] = ndegrees++;
             }
             else 
               degrees[phtable[to]] += adjwgt[j];
           }
        }

        /* Determine which of the ndegrees moves is the best */
        for (jbest=-1, j=0; j<ndegrees; j++) {
           if (jbest == -1 || degrees[j] > degrees[jbest]) 
             jbest = j;
        }
        if (jbest != -1 && degrees[jbest] < id)
          jbest = -1;


        if (jbest != -1) {
          to = ptarget[jbest];
          old = ARATIO(dim, psurf[from], pvol[from]) + ARATIO(dim, psurf[to], pvol[to]);
          new = ARATIO(dim, psurf[from]+id-ed-vsurf[i], pvol[from]-vvol[i]) +
                ARATIO(dim, psurf[to]+ed+id-2.0*degrees[jbest]+vsurf[i], pvol[to]+vvol[i]);
          best = old-new;

          if (best >= 0.0 || degrees[jbest]-id + best > 0.0) {
            where[i] = to;
            INC_DEC(pwgts[to], pwgts[from], vwgt[i]); 
            INC_DEC(pvol[to], pvol[from], vvol[i]);
            psurf[from] = psurf[from] + id - ed - vsurf[i];
            psurf[to]   = psurf[to] + id + ed - 2.0*degrees[jbest] + vsurf[i];
            graph->minratio = ComputeFunction(ctrl->RType, ctrl, graph);
            nmoves++;

            IFSET(ctrl->dbglvl, DBG_MOVEINFO,
              printf("\tMoving %6d from %3d to %3d. Gain: %4.2f. MinRatio: %e [%e]\n",
                     i, from, to, best, graph->minratio, vsurf[i]));
          }
          /* CheckParams(ctrl, graph); */
        }

        for (j=0; j<ndegrees; j++) 
           phtable[ptarget[j]] = -1;
     }

     IFSET(ctrl->dbglvl, DBG_REFINE,
       printf("\t[%6d %6d], Nmoves: %5d, MinRatio: %e\n",
               pwgts[iamin(nparts, pwgts)], pwgts[iamax(nparts, pwgts)], nmoves,
               graph->minratio));

     if (nmoves == 0)
       break;
  }

  graph->nmoves = nmoves;
  IMfree(&perm, &phtable, &degrees, &ptarget, LTERM);
}

/*************************************************************************
* This function performs k-way refinement, whose objective is to directly
* minimize the maximum aspect ratio.
* It also preforms moves that "average" the aspect ratios
* between two domains.
**************************************************************************/
void Random_KWayMinMaxAverageARatioRefine(CtrlType *ctrl, GraphType *graph,
                                          int npasses)
{
  int i, ii, j, dim, nparts, pass, nvtxs, nmoves, ndegrees, pmax; 
  int from, to, jbest;
  idxtype *xadj, *vwgt, *adjncy, *where, *pwgts, *perm, *phtable, *ptarget;
  realtype old, new, best, id, ed, maxar;
  realtype OldToAR, NewToAR, OldFromAR, NewFromAR;
  realtype *vvol, *vsurf, *adjwgt, *pvol, *psurf, *degrees;

  nvtxs     = graph->nvtxs;
  xadj      = graph->xadj;
  vwgt      = graph->vwgt;
  vvol      = graph->vvol;
  vsurf     = graph->vsurf;
  adjncy    = graph->adjncy;
  adjwgt    = graph->adjwgt;

  where  = graph->where;
  pwgts  = graph->pwgts;
  pvol   = graph->pvol;
  psurf  = graph->psurf;
  
  dim      = ctrl->dim;
  nparts   = ctrl->nparts;
  degrees  = realmalloc(nparts, "degrees");
  phtable  = idxsmalloc(nparts, -1, "phtable");
  ptarget  = idxsmalloc(nparts, -1, "ptarget");
  perm     = idxmalloc(nvtxs, "perm");

  /* Determine the domain that has the maximum aspect ratio */
  pmax = 0;
  maxar = ARATIO(dim, psurf[0], pvol[0]);
  for (i=1; i<nparts; i++) {
     new = ARATIO(dim, psurf[i], pvol[i]);
     if (new > maxar) {
       maxar = new;
       pmax = i;
     }
  }

  IFSET(ctrl->dbglvl, DBG_REFINE,
     printf("Partitions: [%3d %3d]-[%3d %3d]. MaxRatio: [%4d, %e], Ratio: %e\n",
             pwgts[iamin(nparts, pwgts)], pwgts[iamax(nparts, pwgts)], 
             ctrl->minsize, ctrl->maxsize, pmax, maxar, graph->minratio));

  RandomPermute(nvtxs, perm, 1);
  for (pass=0; pass<npasses; pass++) {
     RandomPermute(nvtxs, perm, 0);
     RandomPermute(nvtxs, perm, 0);
     for (nmoves=ii=0; ii<nvtxs; ii++) {
        i = perm[ii];
        from = where[i];

        if (pwgts[from] - vwgt[i] < ctrl->minsize)
          continue;

        /* Determine the connectivity of the 'i' vertex */
        for (id=ed=0.0, ndegrees=0, j=xadj[i]; j<xadj[i+1]; j++) {
           to = where[adjncy[j]];

           if (to == from)
             id += adjwgt[j];
           else
             ed += adjwgt[j];

           if (to != from && pwgts[to]+vwgt[i] <= ctrl->maxsize) {
             if (phtable[to] == -1) {
               degrees[ndegrees] = adjwgt[j];
               ptarget[ndegrees] = to;
               phtable[to] = ndegrees++;
             }
             else 
               degrees[phtable[to]] += adjwgt[j];
           }
        }

        /* Determine which of the ndegrees moves is the best */
        for (best=0.01, jbest=-1, j=0; j<ndegrees; j++) {
           to = ptarget[j];
           OldFromAR = ARATIO(dim, psurf[from], pvol[from]);
           OldToAR   = ARATIO(dim, psurf[to], pvol[to]);
           NewFromAR = ARATIO(dim, psurf[from]+id-ed-vsurf[i], pvol[from]-vvol[i]);
           NewToAR   = ARATIO(dim, psurf[to]+ed+id-2.0*degrees[j]+vsurf[i], pvol[to]+vvol[i]);

           /* Check if it increases the max aspect ratio */
           if (NewFromAR > maxar || NewToAR > maxar)
             continue;

           /* If not... */
           /* If move involves partition with max asp ratio, do the move now */
           if (to == pmax || from == pmax) {
             jbest = j;
             break;
           }

           /* Else if partition with max asp ratio is not involved, do the move
              that gives best local gain */
           else {
             old = amax (OldFromAR, OldToAR);
             new = amax (NewFromAR, NewToAR);

             if (old-new > best) {
               best = old-new;
               jbest = j;
             }
           }
        }

        if (jbest != -1) {
          to = ptarget[jbest];

          where[i] = to;
          INC_DEC(pwgts[to], pwgts[from], vwgt[i]); 
          INC_DEC(pvol[to], pvol[from], vvol[i]);
          psurf[from] = psurf[from] + id - ed - vsurf[i];
          psurf[to]   = psurf[to] + id + ed - 2.0*degrees[jbest] + vsurf[i];
     
          /* If we moved from/to the pmax subdomain find the new one ! */
          if (from == pmax || to == pmax) {
            pmax = 0;
            maxar = ARATIO(dim, psurf[0], pvol[0]);
            for (i=1; i<nparts; i++) {
               new = ARATIO(dim, psurf[i], pvol[i]);
               if (new > maxar) {
                 maxar = new;
                 pmax = i;
               }
            }
     
            graph->minratio = maxar;
          }

          nmoves++;

          IFSET(ctrl->dbglvl, DBG_MOVEINFO,
             printf("\tMoving %6d from %3d to %3d. Gain: %4.2f. MinRatio: %e [%e]\n",
                          i, from, to, best, graph->minratio, vsurf[i]));

          /* CheckParams(ctrl, graph); */
        }

        for (j=0; j<ndegrees; j++) 
           phtable[ptarget[j]] = -1;
     }

     IFSET(ctrl->dbglvl, DBG_REFINE, printf("\t[%6d %6d], Nmoves: %5d, MinRatio: %e\n",
                     pwgts[iamin(nparts, pwgts)], pwgts[iamax(nparts, pwgts)], 
                     nmoves, graph->minratio));

     if (nmoves == 0)
       break;
  }

  graph->nmoves = nmoves;
  IFSET(ctrl->dbglvl, DBG_REFINE, printf("FinalMax: %d %e\n", pmax, maxar));

  IMfree(&perm, &phtable, &degrees, &ptarget, LTERM);
}


/*************************************************************************
* This function performs k-way refinement, whose objective is to directly
* minimize the maximum aspect ratio.
**************************************************************************/
void Random_KWayMinMaxARatioRefine(CtrlType *ctrl, GraphType *graph, int npasses)
{
  int i, ii, j, dim, nparts, pass, nvtxs, nmoves, ndegrees, pmax; 
  int from, to, jbest;
  idxtype *xadj, *vwgt, *adjncy, *where, *pwgts, *perm, *phtable, *ptarget;
  realtype new, best, id, ed, maxar;
  realtype NewToAR, NewFromAR;
  realtype *vvol, *vsurf, *adjwgt, *pvol, *psurf, *degrees;

  nvtxs     = graph->nvtxs;
  xadj      = graph->xadj;
  vwgt      = graph->vwgt;
  vvol      = graph->vvol;
  vsurf     = graph->vsurf;
  adjncy    = graph->adjncy;
  adjwgt    = graph->adjwgt;

  where  = graph->where;
  pwgts  = graph->pwgts;
  pvol   = graph->pvol;
  psurf  = graph->psurf;
  
  dim      = ctrl->dim;
  nparts   = ctrl->nparts;
  degrees  = realmalloc(nparts, "degrees");
  phtable  = idxsmalloc(nparts, -1, "phtable");
  ptarget  = idxsmalloc(nparts, -1, "ptarget");
  perm     = idxmalloc(nvtxs, "perm");

  /* Determine the domain that has the maximum aspect ratio */
  pmax = 0;
  maxar = ARATIO(dim, psurf[0], pvol[0]);
  for (i=1; i<nparts; i++) {
     new = ARATIO(dim, psurf[i], pvol[i]);
     if (new > maxar) {
       maxar = new;
       pmax = i;
     }
  }

  IFSET(ctrl->dbglvl, DBG_REFINE,
     printf("Partitions: [%3d %3d]-[%3d %3d]. MaxRatio: [%4d, %e], Ratio: %e\n",
             pwgts[iamin(nparts, pwgts)], pwgts[iamax(nparts, pwgts)], 
             ctrl->minsize, ctrl->maxsize, pmax, maxar, graph->minratio));

  RandomPermute(nvtxs, perm, 1);
  for (pass=0; pass<npasses; pass++) {
     RandomPermute(nvtxs, perm, 0);
     RandomPermute(nvtxs, perm, 0);
     for (nmoves=ii=0; ii<nvtxs; ii++) {
        i = perm[ii];
        from = where[i];

        if (pwgts[from] - vwgt[i] < ctrl->minsize)
          continue;

        /* Determine the connectivity of the 'i' vertex */
        for (id=ed=0.0, ndegrees=0, j=xadj[i]; j<xadj[i+1]; j++) {
           to = where[adjncy[j]];

           if (to == from)
             id += adjwgt[j];
           else
             ed += adjwgt[j];

           if (to != from && pwgts[to]+vwgt[i] <= ctrl->maxsize) {
             if (phtable[to] == -1) {
               degrees[ndegrees] = adjwgt[j];
               ptarget[ndegrees] = to;
               phtable[to] = ndegrees++;
             }
             else 
              degrees[phtable[to]] += adjwgt[j];
           }
        }

        /* Determine which of the ndegrees moves is the best */
        for (jbest=-1, j=0; j<ndegrees; j++) {
           to = ptarget[j];
           NewFromAR = ARATIO(dim, psurf[from]+id-ed-vsurf[i], pvol[from]-vvol[i]);
           NewToAR = ARATIO(dim, psurf[to]+ed+id-2.0*degrees[j]+vsurf[i], pvol[to]+vvol[i]);

           /* Check if move involves partition with max asp ratio and */
           /* if it decreases the max aspect ratio */
           if ( (to == pmax || from == pmax) && NewFromAR <= maxar && NewToAR <= maxar) {
             jbest = j;
             break;
           }
        }

        if (jbest != -1) {
          to = ptarget[jbest];
          where[i] = to;
          INC_DEC(pwgts[to], pwgts[from], vwgt[i]); 
          INC_DEC(pvol[to], pvol[from], vvol[i]);
          psurf[from] = psurf[from] + id - ed - vsurf[i];
          psurf[to]   = psurf[to] + id + ed - 2.0*degrees[jbest] + vsurf[i];
          best = maxar;
          nmoves++;

          /* If we moved to/from to the pmax domain find the new one! */
          if (from == pmax || to==pmax) {
            pmax = 0;
            maxar = ARATIO(dim, psurf[0], pvol[0]);
            for (i=1; i<nparts; i++) {
               new = ARATIO(dim, psurf[i], pvol[i]);
               if (new > maxar) {
                 maxar = new;
                 pmax = i;
               }
            }

            graph->minratio = maxar;
          }

          best -= maxar;
          IFSET(ctrl->dbglvl, DBG_MOVEINFO,
             printf("\tMoving %6d from %3d to %3d. Gain: %4.2f. MinRatio: %e [%e]\n",                    i, from, to, best, graph->minratio, vsurf[i]));
          /* CheckParams(ctrl, graph); */
        }

        for (j=0; j<ndegrees; j++) 
           phtable[ptarget[j]] = -1;
     }

     IFSET(ctrl->dbglvl, DBG_REFINE,
                 printf("\t[%6d %6d], Nmoves: %5d, MinRatio: %e\n",
                 pwgts[iamin(nparts, pwgts)], pwgts[iamax(nparts, pwgts)],
                 nmoves, graph->minratio));

     if (nmoves == 0)
       break;
  }

  graph->nmoves = nmoves;
  IFSET(ctrl->dbglvl, DBG_REFINE, printf("FinalMax: %d %e\n", pmax, maxar));

  IMfree(&perm, &phtable, &degrees, &ptarget, LTERM);
}


/*************************************************************************
* This function performs k-way refinement, whose objective is to directly
* minimize the multiple objective aspect ratio
* In this case RType = REFINE_MINMAXAVAR and RType =  REFINE_WAR
**************************************************************************/
void Random_KWayMultiObjRefine(CtrlType *ctrl, GraphType *graph, int npasses)
{
  int i, ii, j, dim, nparts, pass, nvtxs, nmoves, ndegrees, pmax; 
  int from, to, jbest, jbest1, jbest2;
  idxtype *xadj, *vwgt, *adjncy, *where, *pwgts, *perm, *phtable, *ptarget;
  realtype old, new, best, best1, best2, id, ed, maxar;
  realtype OldToAR, NewToAR, OldFromAR, NewFromAR;
  realtype *vvol, *vsurf, *adjwgt, *pvol, *psurf, *degrees;

  nvtxs     = graph->nvtxs;
  xadj      = graph->xadj;
  vwgt      = graph->vwgt;
  vvol      = graph->vvol;
  vsurf     = graph->vsurf;
  adjncy    = graph->adjncy;
  adjwgt    = graph->adjwgt;

  where  = graph->where;
  pwgts  = graph->pwgts;
  pvol   = graph->pvol;
  psurf  = graph->psurf;
  
  dim      = ctrl->dim;
  nparts   = ctrl->nparts;
  degrees  = realmalloc(nparts, "degrees");
  phtable  = idxsmalloc(nparts, -1, "phtable");
  ptarget  = idxsmalloc(nparts, -1, "ptarget");
  perm     = idxmalloc(nvtxs, "perm");

  /* Determine the domain that has the maximum aspect ratio */
  pmax = 0;
  maxar = ARATIO(dim, psurf[0], pvol[0]);
  for (i=1; i<nparts; i++) {
     new = ARATIO(dim, psurf[i], pvol[i]);
     if (new > maxar) {
       maxar = new;
       pmax = i;
     }
  }

  IFSET(ctrl->dbglvl, DBG_REFINE,
     printf("Partitions: [%3d %3d]-[%3d %3d]. MaxRatio: [%4d, %e], Ratio: %e\n",
             pwgts[iamin(nparts, pwgts)], pwgts[iamax(nparts, pwgts)], 
             ctrl->minsize, ctrl->maxsize, pmax, maxar, graph->minratio));

  RandomPermute(nvtxs, perm, 1);

  for (pass=0; pass<npasses; pass++) {
     RandomPermute(nvtxs, perm, 0);
     RandomPermute(nvtxs, perm, 0);

     for (nmoves=ii=0; ii<nvtxs; ii++) {
        i = perm[ii];
        from = where[i];

        if (pwgts[from] - vwgt[i] < ctrl->minsize)
          continue;

        /* Determine the connectivity of the 'i' vertex */
        for (id=ed=0.0, ndegrees=0, j=xadj[i]; j<xadj[i+1]; j++) {
           to = where[adjncy[j]];

           if (to == from)
             id += adjwgt[j];
           else
             ed += adjwgt[j];

           if (to != from && pwgts[to]+vwgt[i] <= ctrl->maxsize) {
             if (phtable[to] == -1) {
               degrees[ndegrees] = adjwgt[j];
               ptarget[ndegrees] = to;
               phtable[to] = ndegrees++;
             }
             else 
               degrees[phtable[to]] += adjwgt[j];
           }
        }

        /* Determine which of the ndegrees moves is the best */
        for (best1=0.01, best2=0.1, jbest1=-1, jbest2=-1, j=0; j<ndegrees; j++) {
           to = ptarget[j];
           OldFromAR = ARATIO(dim, psurf[from], pvol[from]);
           OldToAR   = ARATIO(dim, psurf[to], pvol[to]);
           NewFromAR = ARATIO(dim, psurf[from]+id-ed-vsurf[i], pvol[from]-vvol[i]);
           NewToAR   = ARATIO(dim, psurf[to]+ed+id-2.0*degrees[j]+vsurf[i], pvol[to]+vvol[i]);

           /* Check first objective min(max) */

           /* Check if it increases the max aspect ratio */
           if (NewFromAR > maxar || NewToAR > maxar)
             continue;

           /* If not... */
           /* If move involves partition with max asp ratio, do the move now */
           if (to == pmax || from == pmax) {
             jbest1 = j;
             break;
           }

           /* Else if partition with max asp ratio is not involved, do the move
              that gives best local gain */
           else {
             old = amax(OldFromAR, OldToAR);
             new = amax(NewFromAR, NewToAR);
             if (old-new > best1) {
               best1  = old-new;
               jbest1 = j;
             }
           }

           /* Check second objective weighted sum */

           old = OldFromAR*pwgts[from] + OldToAR*pwgts[to];
           new = NewFromAR*(pwgts[from]-vwgt[i]) + NewToAR*(pwgts[to]+vwgt[i]);

           if (best2 < old-new) {
             best2 = old-new;
             jbest2 = j;
           }
        }

        IFSET(ctrl->dbglvl, DBG_MOVEINFO,
              printf("\tjbest1=%d, jbest2=%d Gains: %8.6f %8.6f.\n", jbest1,
                     jbest2, best1, best2));

        if (jbest1 != -1)  {
          jbest = jbest1;
          best = best1;
          IFSET(ctrl->dbglvl, DBG_MOVEINFO,
                printf("\t1st OBJECTIVE. Gain: %8.6f\n", best1));
        }  
        else if (jbest2 != -1) {
          jbest = jbest2;
          best = best2;
          IFSET(ctrl->dbglvl, DBG_MOVEINFO,
                printf("\t2nd OBJECTIVE. Gain: %8.6f\n", best2));
        }
        else {
          jbest = -1;
          IFSET(ctrl->dbglvl, DBG_MOVEINFO,
                printf("\tNO OBJECTIVE. Gains: %8.6f %8.6f.\n" , best1, best2));
        }

        if (jbest != -1) {
          to = ptarget[jbest];

          where[i] = to;
          INC_DEC(pwgts[to], pwgts[from], vwgt[i]); 
          INC_DEC(pvol[to], pvol[from], vvol[i]);
          psurf[from] = psurf[from] + id - ed - vsurf[i];
          psurf[to]   = psurf[to] + id + ed - 2.0*degrees[jbest] + vsurf[i];

          /* If we moved from/to the pmax subdomain find the new one! */
          if (from == pmax || to == pmax) {
            pmax = 0;
            maxar = ARATIO(dim, psurf[0], pvol[0]);
            for (i=1; i<nparts; i++) {
               if ((new = ARATIO(dim, psurf[i], pvol[i])) > maxar) {
                 maxar = new;
                 pmax = i;
               }
            }

            graph->minratio = maxar;
          }
          nmoves++;

          IFSET(ctrl->dbglvl, DBG_MOVEINFO,
              printf("\tMoving %6d from %3d to %3d. Gain: %4.2f. MinRatio: %e [%e]\n"                    , i, from, to, best, graph->minratio, vsurf[i]));

          /* CheckParams(ctrl, graph); */
        }

        for (j=0; j<ndegrees; j++) 
          phtable[ptarget[j]] = -1;
     }

     IFSET(ctrl->dbglvl, DBG_REFINE,
           printf("\t[%6d %6d], Nmoves: %5d, MinRatio: %e\n",
           pwgts[iamin(nparts, pwgts)], pwgts[iamax(nparts, pwgts)],
           nmoves, graph->minratio));

     if (nmoves == 0)
       break;
  }

  graph->nmoves = nmoves;
  IFSET(ctrl->dbglvl, DBG_REFINE, printf("FinalMax: %d %e\n", pmax, maxar));

  IMfree(&perm, &phtable, &degrees, &ptarget, LTERM);
}


/*************************************************************************
* This function performs k-way refinement, whose objective is to directly
* minimize the multiple objective aspect ratio
* In this case RType = REFINE_MINMAXAR and RType =  REFINE_WAR
**************************************************************************/
void Random_KWayMultiObjRefine2(CtrlType *ctrl, GraphType *graph, int npasses)
{
  int i, ii, j, dim, nparts, pass, nvtxs, nmoves, ndegrees, pmax; 
  int from, to, jbest, jbest1, jbest2;
  idxtype *xadj, *vwgt, *adjncy, *where, *pwgts, *perm, *phtable, *ptarget;
  realtype old, new, best, best1, best2, id, ed, maxar;
  realtype OldToAR, NewToAR, OldFromAR, NewFromAR;
  realtype *vvol, *vsurf, *adjwgt, *pvol, *psurf, *degrees;

  nvtxs     = graph->nvtxs;
  xadj      = graph->xadj;
  vwgt      = graph->vwgt;
  vvol      = graph->vvol;
  vsurf     = graph->vsurf;
  adjncy    = graph->adjncy;
  adjwgt    = graph->adjwgt;

  where  = graph->where;
  pwgts  = graph->pwgts;
  pvol   = graph->pvol;
  psurf  = graph->psurf;
  
  dim      = ctrl->dim;
  nparts   = ctrl->nparts;
  degrees  = realmalloc(nparts, "degrees");
  phtable  = idxsmalloc(nparts, -1, "phtable");
  ptarget  = idxsmalloc(nparts, -1, "ptarget");
  perm     = idxmalloc(nvtxs, "perm");

  /* Determine the domain that has the maximum aspect ratio */
  pmax = 0;
  maxar = ARATIO(dim, psurf[0], pvol[0]);
  for (i=1; i<nparts; i++) {
     new = ARATIO(dim, psurf[i], pvol[i]);
     if (new > maxar) {
       maxar = new;
       pmax = i;
     }
  }

  IFSET(ctrl->dbglvl, DBG_REFINE,
     printf("Partitions: [%3d %3d]-[%3d %3d]. MaxRatio: [%4d, %e], Ratio: %e\n",
             pwgts[iamin(nparts, pwgts)], pwgts[iamax(nparts, pwgts)], 
             ctrl->minsize, ctrl->maxsize, pmax, maxar, graph->minratio));

  RandomPermute(nvtxs, perm, 1);

  for (pass=0; pass<npasses; pass++) {
     RandomPermute(nvtxs, perm, 0);
     RandomPermute(nvtxs, perm, 0);

     for (nmoves=ii=0; ii<nvtxs; ii++) {
        i = perm[ii];
        from = where[i];

        if (pwgts[from] - vwgt[i] < ctrl->minsize)
          continue;

        /* Determine the connectivity of the 'i' vertex */
        for (id=ed=0.0, ndegrees=0, j=xadj[i]; j<xadj[i+1]; j++) {
           to = where[adjncy[j]];

           if (to == from)
             id += adjwgt[j];
           else
             ed += adjwgt[j];

           if (to != from && pwgts[to]+vwgt[i] <= ctrl->maxsize) {
             if (phtable[to] == -1) {
               degrees[ndegrees] = adjwgt[j];
               ptarget[ndegrees] = to;
               phtable[to] = ndegrees++;
             }
             else 
               degrees[phtable[to]] += adjwgt[j];
           }
        }

        /* Determine which of the ndegrees moves is the best */
        for (best1=0.01, best2=0.1, jbest1=-1, jbest2=-1, j=0; j<ndegrees; j++) {
           to = ptarget[j];
           OldFromAR = ARATIO(dim, psurf[from], pvol[from]);
           OldToAR   = ARATIO(dim, psurf[to], pvol[to]);
           NewFromAR = ARATIO(dim, psurf[from]+id-ed-vsurf[i], pvol[from]-vvol[i]);
           NewToAR   = ARATIO(dim, psurf[to]+ed+id-2.0*degrees[j]+vsurf[i], pvol[to]+vvol[i]);

           /* Check first objective min(max) */

           /* Check if it increases the max aspect ratio */
           if (NewFromAR > maxar || NewToAR > maxar)
             continue;

           /* If not... */
           /* If move involves partition with max asp ratio, do the move now */
           if (to == pmax || from == pmax) {
             jbest1 = j;
             break;
           }

           /* Check second objective weighted sum */

           old = OldFromAR*pwgts[from] + OldToAR*pwgts[to];
           new = NewFromAR*(pwgts[from]-vwgt[i]) + NewToAR*(pwgts[to]+vwgt[i]);

           if (best2 < old-new) {
             best2 = old-new;
             jbest2 = j;
           }
        }

        IFSET(ctrl->dbglvl, DBG_MOVEINFO,
              printf("\tjbest1=%d, jbest2=%d Gains: %8.6f %8.6f.\n", jbest1,
                     jbest2, best1, best2));

        if (jbest1 != -1)  {
          jbest = jbest1;
          best = best1;
          IFSET(ctrl->dbglvl, DBG_MOVEINFO,
                printf("\t1st OBJECTIVE. Gain: %8.6f\n", best1));
        }  
        else if (jbest2 != -1) {
          jbest = jbest2;
          best = best2;
          IFSET(ctrl->dbglvl, DBG_MOVEINFO,
                printf("\t2nd OBJECTIVE. Gain: %8.6f\n", best2));
        }
        else {
          jbest = -1;
          IFSET(ctrl->dbglvl, DBG_MOVEINFO,
                printf("\tNO OBJECTIVE. Gains: %8.6f %8.6f.\n" , best1, best2));
        }

        if (jbest != -1) {
          to = ptarget[jbest];

          where[i] = to;
          INC_DEC(pwgts[to], pwgts[from], vwgt[i]); 
          INC_DEC(pvol[to], pvol[from], vvol[i]);
          psurf[from] = psurf[from] + id - ed - vsurf[i];
          psurf[to]   = psurf[to] + id + ed - 2.0*degrees[jbest] + vsurf[i];

          /* If we moved from/to the pmax subdomain find the new one! */
          if (from == pmax || to == pmax) {
            pmax = 0;
            maxar = ARATIO(dim, psurf[0], pvol[0]);
            for (i=1; i<nparts; i++) {
               if ((new = ARATIO(dim, psurf[i], pvol[i])) > maxar) {
                 maxar = new;
                 pmax = i;
               }
            }

            graph->minratio = maxar;
          }
          nmoves++;

          IFSET(ctrl->dbglvl, DBG_MOVEINFO,
              printf("\tMoving %6d from %3d to %3d. Gain: %4.2f. MinRatio: %e [%e]\n"                    , i, from, to, best, graph->minratio, vsurf[i]));

          /* CheckParams(ctrl, graph); */
        }

        for (j=0; j<ndegrees; j++) 
          phtable[ptarget[j]] = -1;
     }

     IFSET(ctrl->dbglvl, DBG_REFINE,
           printf("\t[%6d %6d], Nmoves: %5d, MinRatio: %e\n",
           pwgts[iamin(nparts, pwgts)], pwgts[iamax(nparts, pwgts)],
           nmoves, graph->minratio));

     if (nmoves == 0)
       break;
  }

  graph->nmoves = nmoves;
  IFSET(ctrl->dbglvl, DBG_REFINE, printf("FinalMax: %d %e\n", pmax, maxar));

  IMfree(&perm, &phtable, &degrees, &ptarget, LTERM);
}


/*************************************************************************
* This function checks the psurf and pvol parameters of the partitioning
**************************************************************************/
void CheckParams(CtrlType *ctrl, GraphType *graph)
{
  int i, j, me, nvtxs, nparts; 
  idxtype *xadj, *vwgt, *adjncy, *where, *pwgts;
  realtype *vvol, *vsurf, *adjwgt, *pvol, *psurf;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  vvol = graph->vvol;
  vsurf = graph->vsurf;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  where = graph->where;

  nparts = ctrl->nparts;
  pwgts = idxsmalloc(nparts, 0, "pwgts");
  pvol = realsmalloc(nparts, 0.0, "pvol");
  psurf = realsmalloc(nparts, 0.0, "psurf");


  for (i=0; i<nvtxs; i++) {
     me = where[i];
     pwgts[me] += vwgt[i];
     pvol[me]  += vvol[i];
     psurf[me] += vsurf[i];

     for (j=xadj[i]; j<xadj[i+1]; j++) {
        if (where[adjncy[j]] != me)
          psurf[me] += adjwgt[j];
     }
  }

  for (i=0; i<nparts; i++) {
     if (pwgts[i] != graph->pwgts[i])
       printf("pwgts: %d %d %d\n", i, pwgts[i], graph->pwgts[i]);

     if (fabs(pvol[i] - graph->pvol[i]) > .01)
       printf("pvol: %d %e %e\n", i, pvol[i], graph->pvol[i]);

     if (fabs(psurf[i] - graph->psurf[i]) > .01) 
       printf("psurf: %d %e %e\n", i, psurf[i], graph->psurf[i]);

  }

  IMfree(&pwgts, &pvol, &psurf, LTERM);
}
