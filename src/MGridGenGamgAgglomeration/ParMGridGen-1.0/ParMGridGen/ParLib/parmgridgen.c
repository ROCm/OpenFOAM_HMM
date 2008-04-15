/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * parmgridgen.c
 *
 * This is the entry point of ParMETIS_AspectRatio
 *
 * George Irene
 */

#include "parmgridgen.h"

/***********************************************************************************
* This function is the entry point of the parallel coarse grid construction.
* This function assumes nothing about the graph distribution.
* It is the general case.
************************************************************************************/
void ParMGridGen(idxtype *vtxdist, idxtype *xadj, realtype *vvol, realtype *vsurf,
                 idxtype *adjncy, realtype *adjwgt, int *nparts, int minsize,
                 int maxsize, int *options, idxtype *part, MPI_Comm *comm)
{
  int i,j;
  int npes, mype;
  int nmoves, nsteps, firstfvtx;
  int orgnvtxs, nvtxs, nedges;
  int moptions[10];
  idxtype **p_glblvtxid, *glblvtxid, *fvtxdist, *perm;
  idxtype *myvtxdist, *myxadj, *myadjncy, *mypart, *where;
  idxtype **p_xadj, **p_adjncy, **p_part;
  idxtype *orxadj, *oradjncy;
  realtype *myvvol, *myvsurf, *myadjwgt;
  realtype **p_vvol, **p_vsurf, **p_adjwgt; 
  realtype *oradjwgt, *orvsurf;
  MGridCtrlType ctrl;
  MGridGraphType *graph, *mgraph;
  MGridWorkSpaceType wspace;
  double tmr1;

  MPI_Comm_size(*comm, &npes);
  MPI_Comm_rank(*comm, &mype);

  IFSET(options[OPTION_DBGLVL], DBG_STEPS, printf("---------------- CYCLE 0 ----------------\n"));
  
  orgnvtxs = nvtxs = vtxdist[mype+1]-vtxdist[mype];
  nedges = xadj[nvtxs];

  /* Create my own copy of the graph, so that the original can be returned intact */
  myvtxdist= idxmalloc(npes+1, "ParMGridGen: myvtxdist");
  myxadj   = idxmalloc(nvtxs+1, "ParMGridGen: myxadj");
  myvvol   = realmalloc(nvtxs, "ParMGridGen: myvvol");
  myvsurf  = realmalloc(nvtxs, "ParMGridGen: myvsurf");
  myadjncy = idxmalloc(nedges, "ParMGridGen: myadjncy");
  myadjwgt = realmalloc(nedges, "ParMGridGen: myadjwgt");
  mypart = idxmalloc(nvtxs, "ParMGridGen: mypart");

  idxcopy(npes+1, vtxdist, myvtxdist);
  idxcopy(nvtxs+1, xadj, myxadj);
  realcopy(nvtxs, vvol, myvvol);
  realcopy(nvtxs, vsurf, myvsurf);
  idxcopy(nedges, adjncy, myadjncy);
  realcopy(nedges, adjwgt, myadjwgt);

  p_xadj = &myxadj;
  p_vvol = &myvvol;
  p_vsurf = &myvsurf;
  p_adjncy = &myadjncy;
  p_adjwgt = &myadjwgt;
  p_part = &mypart;
  
/*
  if (*numflag == 1) 
    ChangeNumbering(myvtxdist, myxadj, myadjncy, mypart, npes, mype, 1);
*/

  /* Keep a copy of the initial graph. */
  /* I will need those twice in the next step when it will have been ruined */
  orxadj   = idxmalloc(nvtxs+1, "ParMGridGen: orxadj");
  oradjncy = idxmalloc(nedges, "ParMGridGen: oradjncy");
  oradjwgt = realmalloc(nedges, "ParMGridGen: oradjwgt");
  orvsurf  = realmalloc(nvtxs, "ParMGridGen: orvsurf");
  
  idxcopy(nvtxs+1, myxadj, orxadj);
  idxcopy(nedges, myadjncy, oradjncy);
  realcopy(nedges, myadjwgt, oradjwgt);
  realcopy(nvtxs, myvsurf, orvsurf);

  glblvtxid = idxmalloc(nvtxs, "ParMGridGen: glblvtxid");
  for (j=0, i=myvtxdist[mype]; i<myvtxdist[mype+1]; i++)
     glblvtxid[j++]=i;
  p_glblvtxid = &glblvtxid;

  /* Create purely local graph by removing edges and correcting vsurf */
  RemoveCorrectSurfaceEdges(myvtxdist, myxadj, myadjncy, myadjwgt, orxadj,
                            oradjncy, oradjwgt, orvsurf, *comm);

  /* Call MGridGen to do the partitioning */
  moptions[OPTION_CTYPE]=options[OPTION_CTYPE];
  moptions[OPTION_RTYPE]=options[OPTION_RTYPE];
  if ( options[OPTION_DBGLVL] == DBG_STEPS )
    moptions[OPTION_DBGLVL]=DBG_TRACK;
  else if ( options[OPTION_DBGLVL] == DBG_TRACK )
    moptions[OPTION_DBGLVL]=0;
  else
    moptions[OPTION_DBGLVL]=options[OPTION_DBGLVL];
  moptions[OPTION_DIM]=options[OPTION_DIM];
/*
  cleartimer(tmr1);
  starttimer(tmr1);
*/
  MGridGen(nvtxs, orxadj, myvvol, orvsurf, oradjncy, oradjwgt, minsize, maxsize,
           moptions, &nmoves, nparts, mypart);
/*
  stoptimer(tmr1);
  printf("MGridGen Time: %lf\n", gettimer(tmr1));
*/

  IMfree(&orxadj, &oradjncy, &oradjwgt, &orvsurf, LTERM);

  IFSET(options[OPTION_DBGLVL], DBG_STEPS, printf("%d END OF CYCLE 0 : MGridGen nparts=%d\n",mype, *nparts));
  /* Move graph around and refine for #nsteps steps */
/*
  cleartimer(tmr1);
  starttimer(tmr1);
*/
  nsteps = 3;
  MoveRefine(myvtxdist, p_xadj, p_vvol, p_vsurf, p_adjncy, p_adjwgt, p_glblvtxid,
             nparts, minsize, maxsize, &nmoves, &nsteps, options, p_part, comm);
/*
  stoptimer(tmr1);
  printf("MGridGenRefine Time: %lf\n", gettimer(tmr1));
*/

  glblvtxid = *p_glblvtxid;
  mypart = *p_part;

  /* Communicate number of parts found */
  nvtxs = myvtxdist[mype+1]-myvtxdist[mype];
  fvtxdist = idxmalloc(npes+1, "ParMGridGen: fvtxdist");
  MPI_Allgather((void *)nparts, 1, MPI_INT, (void *)fvtxdist, 1, MPI_INT, *comm);

  MAKECSR(i, npes, fvtxdist);
  firstfvtx = fvtxdist[mype];
  IMfree(&fvtxdist, LTERM);

  /* Convert the part array into global numbering */
  for (i=0; i<nvtxs; i++)
     mypart[i] +=firstfvtx;

  IFSET(options[OPTION_DBGLVL], DBG_TRACK,
               PrintAspectRatioStats(myvtxdist, *p_xadj, *p_vvol, *p_vsurf, *p_adjncy,
                        *p_adjwgt, nparts, minsize, maxsize, options, *p_part, *comm));

  /* Communicate fused element information to correspond to the original graph */
  where = idxsmalloc(nvtxs, -1, "ParMGridGen: sendind");

  for (i=0; i<npes; i++)
     for (j=0; j<nvtxs; j++)
        if (glblvtxid[j] >= vtxdist[i] && glblvtxid[j] < vtxdist[i+1])
          where[j] = i;

  SetUpMGridCtrl(&ctrl, minsize, maxsize, options, *comm);
  ctrl.nparts = npes;

  for (i=0; i<nvtxs; i++)
    ASSERTP(&ctrl, where[i] >= 0 && where[i] < npes, (&ctrl, "%d %d %d\n", i, where[i], npes) );

  graph = SetUpMGridGraph(&ctrl, myvtxdist, *p_xadj, *p_vvol, *p_vsurf, *p_adjncy, *p_adjwgt);
  graph->where = where;
  graph->fusedinfo = mypart;
  graph->glblvtxid = glblvtxid;

  PreAllocateMGridMemory(&ctrl, graph, &wspace);
  
  MGridSetUp(&ctrl, graph, &wspace);

  mgraph = MoveMGridGraph(&ctrl, graph, &wspace);

  nvtxs = mgraph->nvtxs;

  IMfree(&graph->lperm, &graph->peind, &graph->recvptr, &graph->recvind, &graph->sendptr,
         &graph->sendind, &graph->pexadj, &graph->peadjncy, &graph->peadjloc, &graph->imap,
         &graph->vwgt, &graph->adjwgtsum, &graph->where, &graph, LTERM);
  FreeMGridWSpace(&wspace);

  /* Sort so as to go back to original order of input */
  ASSERT(&ctrl, nvtxs == orgnvtxs);
  perm = idxmalloc(nvtxs, "ParMGridGen: perm");

  for (i=0; i<nvtxs; i++)
     perm[mgraph->glblvtxid[i] - vtxdist[mype]] = i;

  for (i=0; i<nvtxs; i++)
     part[i] = mgraph->fusedinfo[perm[i]];
  
  IMfree(&mgraph->vtxdist, &mgraph->xadj, &mgraph->vwgt, &mgraph->vvol, &mgraph->vsurf,
         &mgraph->adjwgtsum, &mgraph->adjncy, &mgraph->adjwgt, &mgraph->fusedinfo,
         &mgraph->glblvtxid, &mgraph,
         p_xadj, p_vvol, p_vsurf, p_adjncy, p_adjwgt, p_part, p_glblvtxid,
         &perm, &myvtxdist, LTERM);
}


/***************************************************************************/
void MoveRefine(idxtype *vtxdist, idxtype **p_xadj, realtype **p_vvol,
                realtype **p_vsurf, idxtype **p_adjncy, realtype **p_adjwgt,
                idxtype **p_glblvtxid, int *nparts, int minsize, int maxsize,
                int *nmoves, int *nsteps, int *options, idxtype **p_part,
                MPI_Comm *comm)
{
  int i, step, count, wgtflag, numflag;
  int firstfvtx, nvtxs, nedges;
  int npes, mype;
  int moptions[10];
  idxtype *mpart, *fusedinfo, *fvtxdist, *glblvtxid;
  idxtype  *xadj, *adjncy, *vwgt, *part;
  idxtype *orxadj, *oradjncy, *oradjncy2;
  realtype *vvol, *adjwgt, *vsurf;
  realtype *oradjwgt, *orvsurf;
  realtype lmeasure1, lmeasure2, gmeasure1, gmeasure2;
  MGridCtrlType ctrl;
  MGridWorkSpaceType wspace;
  MGridGraphType *graph = 0, *mgraph = 0;
  double tmr1;

  MPI_Comm_size(*comm, &npes);
  MPI_Comm_rank(*comm, &mype);

  xadj = *p_xadj;
  vvol = *p_vvol;
  vsurf = *p_vsurf;
  adjncy = *p_adjncy;
  adjwgt = *p_adjwgt;
  part = *p_part;
  glblvtxid = *p_glblvtxid;

  /* Keep initial graph */

  nvtxs = vtxdist[mype+1]-vtxdist[mype];
  nedges = xadj[nvtxs];
  oradjncy2 = idxmalloc(nedges, "CorrectSurfaceEdges: oradjncy2");
  
  idxcopy(nedges, adjncy, oradjncy2); 

  vwgt = idxsmalloc(nvtxs, 1, "MoveRefine: vwgt");

  for (step=0; step<*nsteps; step++) {

     IFSET(options[OPTION_DBGLVL], DBG_STEPS, printf("---------------- CYCLE %d ----------------\n",step+1)); 
     /* Communicate number of parts found in the previous call of MGridgen */ 
     fvtxdist = idxmalloc(npes+1, "ParMETIS_AspectRatio: fvtxdist");
     MPI_Allgather((void *)nparts, 1, MPI_INT, (void *)fvtxdist, 1, MPI_INT, *comm);
   
     MAKECSR(i, npes, fvtxdist);
     firstfvtx = fvtxdist[mype];
   
     /* Convert the part array into global numbering */
     for (i=0; i<nvtxs; i++)
        part[i] +=firstfvtx;
   
     fusedinfo =  idxmalloc(nvtxs, "MoveRefine: fusedinfo"); 
     idxcopy(nvtxs, part, fusedinfo);
   
     /***** Create Fused Element Graph and get new partition *****/
     /* on exit : adjncy array is transformed to local numbering */
     moptions[0]=1;
     wgtflag = 1;
     numflag = 0;
/*  
     cleartimer(tmr1);
     starttimer(tmr1);
*/
     ParMETIS_FusedElementGraph(vtxdist, xadj, vvol, vsurf, adjncy, vwgt, adjwgt,
                       &wgtflag, &numflag, nparts, options, part, comm);
/*
     stoptimer(tmr1);
     printf("ParMETIS_FusedElementGraph Time = %lf\n", gettimer(tmr1)/nsteps);
*/
/*
     count = 0;
     printf("=========================\n");
     for (i=0; i<nvtxs; i++)
        if ( mype != part[i])
          count++;
     lmeasure1 = ((double)(count))/(*nparts);
     lmeasure2 = ((double)(count))/nvtxs;
     printf("count=%d nvtxs=%d nparts=%d count/npart=%f count/nvtxs=%f\n",count,nvtxs, *nparts, lmeasure1, lmeasure2);
     printf("=========================\n");
*/
   
     /* Move graph around based on the results of the fused graph partition */
/*
     cleartimer(tmr1);
     starttimer(tmr1);
*/
     SetUpMGridCtrl(&ctrl, minsize, maxsize, options, *comm); 
/*
     stoptimer(tmr1);
     printf("SetUpMGridCtrl Time = %lf\n",gettimer(tmr1)/nsteps);
*/    
     ctrl.nparts=npes;
   
     IMfree(&graph, &mgraph, LTERM);
/*
     cleartimer(tmr1);
     starttimer(tmr1);
*/
     graph = SetUpMGridGraph(&ctrl, vtxdist, xadj, vvol, vsurf, oradjncy2, adjwgt);
/*
     stoptimer(tmr1);
     printf("SetUpMGridGraph Time = %lf\n",gettimer(tmr1)/nsteps);
*/
     graph->where = part;
     graph->fusedinfo = fusedinfo;
     graph->glblvtxid = glblvtxid;
   
     PreAllocateMGridMemory(&ctrl, graph, &wspace);
   
     MGridSetUp(&ctrl, graph, &wspace);

/*
     cleartimer(tmr1);
     starttimer(tmr1);
*/
     mgraph = MoveMGridGraph(&ctrl, graph, &wspace);
/*
     stoptimer(tmr1);
     printf("MoveMGridGraph Time = %lf\n",gettimer(tmr1)/nsteps);
*/

     IMfree(&xadj, &vvol, &vsurf, &adjwgt, &adjncy, &oradjncy2, &part, &vwgt, 
            &fvtxdist, &fusedinfo, &glblvtxid, &graph->vwgt, &graph->adjwgtsum,
            &graph->lperm, &graph->peind, &graph->recvptr, &graph->recvind,
            &graph->sendind, &graph->sendptr, &graph->pexadj, &graph->peadjncy,
            &graph->peadjloc, &graph->imap, LTERM);

     nvtxs = mgraph->nvtxs;
     nedges = mgraph->xadj[nvtxs];
     orxadj = idxmalloc(nvtxs+1, "CorrectSurfaceEdges: orxadj");
     oradjncy = idxmalloc(nedges, "CorrectSurfaceEdges: oradjncy");
     oradjncy2 = idxmalloc(nedges, "CorrectSurfaceEdges: oradjncy"); 
     oradjwgt = realmalloc(nedges, "CorrectSurfaceEdges: oradjwgt");
     orvsurf = realmalloc(nvtxs, "CorrectSurfaceEdges: orvsurf");
  
     idxcopy(nvtxs+1, mgraph->xadj, orxadj);
     idxcopy(nedges, mgraph->adjncy, oradjncy);
     idxcopy(nedges, mgraph->adjncy, oradjncy2); 
     realcopy(nedges, mgraph->adjwgt, oradjwgt);
     realcopy(nvtxs, mgraph->vsurf, orvsurf);

     mpart = idxmalloc(nvtxs, "MoveRefine: mpart"); 

     /* Create purely local graph by removing edges and correcting vsurf */
/*
     cleartimer(tmr1);
     starttimer(tmr1);
*/
     RemoveCorrectSurfaceEdges(mgraph->vtxdist, orxadj, oradjncy, oradjwgt, 
              mgraph->xadj, mgraph->adjncy, mgraph->adjwgt, mgraph->vsurf, *comm); 
/*
     stoptimer(tmr1); 
     printf("RemoveCorrectSurfaceEdges Time = %lf\n",gettimer(tmr1)/nsteps);
*/
     /* One more step of refinement */
     moptions[OPTION_CTYPE]=options[OPTION_CTYPE];
     moptions[OPTION_RTYPE]=options[OPTION_RTYPE];
     if ( options[OPTION_DBGLVL] == DBG_STEPS )
       moptions[OPTION_DBGLVL]=DBG_TRACK;
     else if ( options[OPTION_DBGLVL] == DBG_TRACK )
       moptions[OPTION_DBGLVL]=0;
     else
       moptions[OPTION_DBGLVL]=options[OPTION_DBGLVL];
     moptions[OPTION_DIM]=options[OPTION_DIM];

/*
     cleartimer(tmr1);
     starttimer(tmr1);
*/
     MGridGenRefine(nvtxs, mgraph->xadj, mgraph->vvol, mgraph->vsurf, mgraph->adjncy, 
                    mgraph->fusedinfo, mgraph->adjwgt, minsize, maxsize, moptions,
                    nmoves, nparts, mpart);

/*
     stoptimer(tmr1);
     printf("MGridGenRefineSerial Time = %lf\n",gettimer(tmr1)/nsteps);
*/

     nedges = mgraph->xadj[nvtxs];

     idxcopy(npes+1, mgraph->vtxdist, vtxdist);
     xadj = orxadj;
     adjncy = oradjncy;
     vvol = mgraph->vvol;
     vsurf = orvsurf;
     vwgt = mgraph->vwgt;
     adjwgt = oradjwgt;
     part = mpart;
     glblvtxid = mgraph->glblvtxid;

/*
     gmeasure1 = MGridGlobalSESumReal(&ctrl, lmeasure1)/npes;
     gmeasure2 = MGridGlobalSESumReal(&ctrl, lmeasure2)/npes;
     printf("count measure1=%f measure2=%f\n",gmeasure1, gmeasure2);
*/

     IMfree(&mgraph->vtxdist, &mgraph->xadj, &mgraph->adjncy, &mgraph->vsurf,
                 &mgraph->adjwgt, &mgraph->adjwgtsum, &mgraph->fusedinfo, LTERM); 
     FreeMGridWSpace(&wspace);
     FreeMGridCtrl(&ctrl);
/*
     if (gmeasure1 < LALALAALA)
       break;
*/
  }

  *p_xadj = xadj;
  *p_vvol = vvol;
  *p_vsurf = vsurf;
  *p_adjncy = adjncy;
  *p_adjwgt = adjwgt;
  *p_part = part;
  *p_glblvtxid = glblvtxid;

  IMfree(&vwgt, &oradjncy2, &graph, &mgraph, LTERM); 
}


/******************************************************************************
* This function takes a graph and its partition vector and then
* computes and prints statistics info regarding the partition found
*******************************************************************************/
void PrintAspectRatioStats(idxtype *vtxdist, idxtype *xadj, realtype *vvol, realtype *vsurf,
                   idxtype *adjncy, realtype *adjwgt, int *nparts, int minsize, int maxsize,
                   int *options, idxtype *part, MPI_Comm comm)
{
  int i, j, k, nvtxs, nedges, dim;
  int firstpart;
  idxtype *myvtxdist, *myxadj, *myadjncy;
  idxtype *lpwgts, *counts;
  idxtype *where, *swhere, *rwhere;
  RInfoType *rinfo, *myrinfo;
  EdgeType *edegrees;
  MGridCtrlType ctrl;
  MGridWorkSpaceType wspace;
  MGridGraphType *graph;
  int me, other;
  int npes, mype;
  realtype *myvvol, *myvsurf, *myadjwgt;
  realtype *lpvol, *lpsurf;
  realtype lmin, lmax, lsum, lwsum, ratio, lsurf;
  realtype gmin, gmax, gsum, gwsum, gsurf;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  dim = options[OPTION_DIM];

  /*IFSET(options[OPTION_DBGLVL], DBG_IRENE, printf("%d PrintAspectRatioStats\n",mype)); */

  nvtxs = vtxdist[mype+1]-vtxdist[mype];
  nedges = xadj[nvtxs];

  /* Create my own copy of the graph, so that the original can be returned intact */
  myvtxdist= idxmalloc(npes+1, "PrintAspectRatioStats: myvtxdist");
  myxadj   = idxmalloc(nvtxs+1, "PrintAspectRatioStats: myxadj");
  myvvol   = realmalloc(nvtxs, "PrintAspectRatioStats: myvvol");
  myvsurf  = realmalloc(nvtxs, "PrintAspectRatioStats: myvsurf");
  myadjncy = idxmalloc(nedges, "PrintAspectRatioStats: myadjncy");
  myadjwgt = realmalloc(nedges, "PrintAspectRatioStats: myadjwgt");

  idxcopy(npes+1, vtxdist, myvtxdist);
  idxcopy(nvtxs+1, xadj, myxadj);
  realcopy(nvtxs, vvol, myvvol);
  realcopy(nvtxs, vsurf, myvsurf);
  idxcopy(nedges, adjncy, myadjncy);
  realcopy(nedges, adjwgt, myadjwgt);

  SetUpMGridCtrl(&ctrl, minsize, maxsize, options, comm);
  ctrl.nparts=*nparts;

  graph = SetUpMGridGraph(&ctrl, myvtxdist, myxadj, myvvol, myvsurf, myadjncy, myadjwgt);

  PreAllocateMGridMemory(&ctrl, graph, &wspace);

  MGridSetUp(&ctrl, graph, &wspace);

  graph->where = idxmalloc(nvtxs+graph->nrecv, "PrintAspectRatioStats: graph->where");
  i = idxamin(nvtxs, part);
  firstpart = part[i];
  for (i=0; i<nvtxs; i++)
     graph->where[i] = part[i] - firstpart;

  counts = idxsmalloc(maxsize+1, 0, "PrintAspectRatioStats: counts");
  where = graph->where;
  rinfo = graph->rinfo = (RInfoType *)IMmalloc(sizeof(RInfoType)*nvtxs, "PrintAspectRatioStats:rinfo");
  lpwgts = graph->lpwgts = idxsmalloc(ctrl.nparts, 0, "PrintAspectRatioStats: lpwgts");
  graph->gpwgts = idxmalloc(ctrl.nparts, "PrintAspectRatioStats: gpwgts");
  lpvol =  graph->lpvol = realsmalloc(ctrl.nparts, 0.0, "PrintAspectRatioStats: lpvol");
  graph->gpvol = realmalloc(ctrl.nparts, "PrintAspectRatioStats: gpvol");
  lpsurf =  graph->lpsurf = realsmalloc(ctrl.nparts, 0.0, "PrintAspectRatioStats: lpsurf");
  graph->gpsurf = realmalloc(ctrl.nparts, "PrintAspectRatioStats: gpsurf");

  /*------------------------------------------------------------
  / Send/Receive the where information of interface vertices
  /------------------------------------------------------------*/
  swhere = wspace.indices;
  rwhere = where + nvtxs;

  MGridCommInterfaceData(&ctrl, graph, where, swhere, rwhere);

  ASSERT(&ctrl, wspace.nlarge >= myxadj[nvtxs]);

  /*------------------------------------------------------------
  / Compute now the id/ed degrees
  /------------------------------------------------------------*/
  graph->lmincut = 0;
  for (i=0; i<nvtxs; i++) {
    me = where[i];
    myrinfo = rinfo+i;

    lpwgts[me] += graph->vwgt[i];
    lpvol[me] += graph->vvol[i];
    lpsurf[me] += graph->vsurf[i];

    myrinfo->degrees = wspace.degrees + myxadj[i];
    myrinfo->ndegrees = 0;
    myrinfo->id = myrinfo->ed = 0.0;

    for (j=myxadj[i]; j<myxadj[i+1]; j++) {
      if (me == where[myadjncy[j]])
/*        myrinfo->id += myadjwgt[j]; */
        myrinfo->id = myrinfo->id + 1.0;
      else {
/*        myrinfo->ed += myadjwgt[j]; */
        myrinfo->ed = myrinfo->ed + 1.0;
        lpsurf[me] += myadjwgt[j];
      }
    }

    if (myrinfo->ed > 0.0) {  /* Time to do some serious work */
      graph->lmincut += myrinfo->ed;
      edegrees = myrinfo->degrees;

      for (j=myxadj[i]; j<myxadj[i+1]; j++) {
        other = where[myadjncy[j]];
        if (me != other) {
          for (k=0; k<myrinfo->ndegrees; k++) {
            if (edegrees[k].edge == other) {
              edegrees[k].ewgt += myadjwgt[j];
              break;
            }
          }
          if (k == myrinfo->ndegrees) {
            edegrees[k].edge = other;
            edegrees[k].ewgt = myadjwgt[j];
            myrinfo->ndegrees++;
          }
          ASSERT(&ctrl, myrinfo->ndegrees <= myxadj[i+1]-myxadj[i]);
        }
      }
    }
  }

  graph->lmincut = graph->lmincut/2;

  lmin = lmax = lsum = ARATIO(dim, lpsurf[0], lpvol[0]);
  lwsum = 1.0*lpwgts[0]*ARATIO(dim, lpsurf[0], lpvol[0]);
  lsurf = lpsurf[0];
  counts[lpwgts[0]]++;
  graph->lminratio = ARATIO(dim, lpsurf[0], lpvol[0]);
  for (i=1; i<ctrl.nparts; i++) {
     ratio = ARATIO(dim, lpsurf[i], lpvol[i]);
     graph->lminratio += ARATIO(dim, lpsurf[i], lpvol[i]);
     lsum += ratio;
     lwsum += 1.0*lpwgts[i]*ratio;
     lsurf += lpsurf[i];
     if (lmin > ratio)
       lmin = ratio;
     if (lmax < ratio)
       lmax = ratio;
     counts[lpwgts[i]]++;
  }

  printf("\n-------------------------- LOCAL --------------------------\n");
  printf("Npoints: %d, Coarsening Factor: %f\n", ctrl.nparts, 1.0*graph->nvtxs/(1.0*ctrl.nparts));
  printf("Aspect Ratios: Min : %e, Max : %e\n", lmin, lmax);
  printf("Aspect Ratios: Sum : %e, Wsum: %e\n", lsum, lwsum);
  printf("Aspect Ratios: Surf: %e, Avg : %e\n", lsurf, lsum/(1.0*ctrl.nparts));
  printf("Graph lmincut : %e\n", graph->lmincut);
  printf("Cell size: min=%d, max=%d\n", ctrl.minsize, ctrl.maxsize);
  for (i=1; i<=maxsize; i++)
     if (counts[i] != 0)
       printf("[%2d %4d] ", i, counts[i]);
  printf("\n");

  /* Finally, sum-up the partition params */
  MPI_Allreduce((void *)&lmin, (void *)&gmin, 1, REAL_DATATYPE, MPI_MIN, comm);
  MPI_Allreduce((void *)&lmax, (void *)&gmax, 1, REAL_DATATYPE, MPI_MAX, comm);
  MPI_Allreduce((void *)&lsum, (void *)&gsum, 1, REAL_DATATYPE, MPI_SUM, comm);
  MPI_Allreduce((void *)&lwsum, (void *)&gwsum, 1, REAL_DATATYPE, MPI_SUM, comm);
  MPI_Allreduce((void *)&lsurf, (void *)&gsurf, 1, REAL_DATATYPE, MPI_SUM, comm);
  graph->mincut = MGridGlobalSESumReal(&ctrl, graph->lmincut);
  i = MGridGlobalSESum(&ctrl, ctrl.nparts);

  printf("-------------------------- GLOBAL --------------------------\n");
  printf("Npoints: %d, Coarsening Factor: %f\n", i, 1.0*graph->gnvtxs/(1.0*i));
  printf("Aspect Ratios: Min : %e, Max : %e\n", gmin, gmax);
  printf("Aspect Ratios: Sum : %e, Wsum: %e\n", gsum, gwsum);
  printf("Aspect Ratios: Surf: %e, Avg : %e\n", gsurf, gsum/(1.0*i));
  printf("Graph mincut : %e\n", graph->mincut);

  printf("\n------------------------------ END ------------------------------\n");

  FreeMGridWSpace(&wspace);
  IMfree(&graph->vtxdist, &graph->xadj, &graph->adjncy, &graph->adjwgt, &graph->vvol,
         &graph->vsurf, &graph->vwgt, &graph->adjwgtsum, &graph->lpvol, &graph->gpvol,
         &graph->lpsurf, &graph->gpsurf, &graph->lperm, &graph->where, &graph->rinfo,
         &graph->lpwgts, &graph->gpwgts, &graph->peind, &graph->sendptr,
         &graph->sendind, &graph->recvptr, &graph->recvind, &graph->imap,
         &graph->pexadj, &graph->peadjncy, &graph->peadjloc, &counts, &graph, LTERM);
}
