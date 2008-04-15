/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mesh2dual.c
 *
 * This file reads in the element node connectivity array of a mesh and writes
 * out its dual in the format suitable for Metis.
 *
 * Started 9/29/97
 * George
 *
 * $Id: mesh2dual.c,v 1.2 2002/08/10 06:02:53 karypis Exp $
 *
 */

#include <metisbin.h>



/*************************************************************************
* Let the game begin
**************************************************************************/
int main(int argc, char *argv[])
{
  idxtype i, j, ne, nn, etype, mtype, cnt, numflag=0;
  idxtype *elmnts, *xadj, *adjncy, *metype;
  idxtype *conmat, *elms, *weights; 
  double IOTmr, DUALTmr;
  char fileout[256], etypestr[5][5] = {"TRI", "TET", "HEX", "QUAD", "LINE"};

  if (argc <2) {
    mprintf("Usage: %s <meshfile> [confile]\n",argv[0]);
    exit(0);
  }


  mtype=MeshType(argv[1]);
  ne=MixedElements(argv[1]);
  metype = idxmalloc(ne, "main: metype");
  weights = idxmalloc(ne, "main: weights");
  
  if (mtype==1 || mtype==3){

  gk_clearcputimer(IOTmr);
  gk_clearcputimer(DUALTmr);

  gk_startcputimer(IOTmr);
  if (mtype==1)
     elmnts = ReadMesh(argv[1], &ne, &nn, &etype);
  else
     elmnts = ReadMeshWgt(argv[1], &ne, &nn, &etype, weights);
  gk_stopcputimer(IOTmr);

  mprintf("**********************************************************************\n");
  mprintf("%s", METISTITLE);
  mprintf("Mesh Information ----------------------------------------------------\n");
  mprintf("  Name: %s, #Elements: %D, #Nodes: %D, Etype: %s\n\n", argv[1], ne, nn, etypestr[etype-1]);
  mprintf("Forming Dual Graph... -----------------------------------------------\n");

  xadj = idxmalloc(ne+1, "main: xadj");
  elms = idxsmalloc(ne+1, 0, "main: elms");


  gk_startcputimer(DUALTmr);
  cnt=METIS_MeshToDualCount(&ne, &nn, elmnts, elms, &etype, &numflag);
  adjncy = idxmalloc(cnt+1, "main: adjncy");
  METIS_MeshToDual(&ne, &nn, elmnts, elms, &etype, &numflag, xadj, adjncy);
  gk_stopcputimer(DUALTmr);

  mprintf("  Dual Information: #Vertices: %D, #Edges: %D\n", ne, xadj[ne]/2);

  msprintf(fileout, "%s.dgraph", argv[1]);
  gk_startcputimer(IOTmr);
  if (mtype==1)
     WriteGraph(fileout, ne, xadj, adjncy);
  else
     WriteWgtGraph(fileout, ne, xadj, adjncy, weights);
     
  gk_stopcputimer(IOTmr);


  mprintf("\nTiming Information --------------------------------------------------\n");
  mprintf("  I/O:          \t\t %7.3f\n", gk_getcputimer(IOTmr));
  mprintf("  Dual Creation:\t\t %7.3f\n", gk_getcputimer(DUALTmr));
  mprintf("**********************************************************************\n");

  }

  else {

  
  gk_clearcputimer(IOTmr);
  gk_clearcputimer(DUALTmr);

  gk_startcputimer(IOTmr);

  
  if(mtype==0)
     elmnts = ReadMixedMesh(argv[1], &ne, &nn, metype);
  else
     elmnts = ReadMixedMeshWgt(argv[1], &ne, &nn, metype, weights);

  if (argc==3)  
  conmat = ReadMgcnums(argv[2]);
  gk_stopcputimer(IOTmr);

   

  mprintf("**********************************************************************\n");
  mprintf("%s", METISTITLE);
  mprintf("Mesh Information ----------------------------------------------------\n");
  mprintf("  Name: %s, #Elements: %D, #Nodes: %D, Etype: %s\n\n", argv[1], ne, nn, "Mixed");
  mprintf("Forming Dual Graph... ----------------------------------------------\n");

  xadj = idxmalloc(ne+1, "main: xadj");
  elms = idxsmalloc(ne+1, 0, "main: elms");

  
  gk_startcputimer(DUALTmr);

  if (argc==3){
  cnt=METIS_MixedMeshToDualCount(&ne, &nn, elmnts, elms, metype, &numflag, 
conmat, 1);
  adjncy = idxmalloc(cnt+1, "main: adjncy");
  METIS_MixedMeshToDual(&ne, &nn, elmnts, elms, metype, &numflag, xadj, adjncy, 
 conmat, 1);
  } 
  else{
  cnt=METIS_MixedMeshToDualCount(&ne, &nn, elmnts, elms, metype, &numflag, 
conmat, 0);
  adjncy = idxmalloc(cnt+1, "main: adjncy");
  METIS_MixedMeshToDual(&ne, &nn, elmnts, elms, metype, &numflag, xadj, adjncy, conmat, 0);
  } 
  gk_stopcputimer(DUALTmr);

  mprintf("  Dual Information: #Vertices: %D, #Edges: %D\n", ne, xadj[ne]/2);

  msprintf(fileout, "%s.dgraph", argv[1]);
  gk_startcputimer(IOTmr);

  if (mtype==0)
     WriteGraph(fileout, ne, xadj, adjncy);
  else
     WriteWgtGraph(fileout, ne, xadj, adjncy, weights);
  gk_stopcputimer(IOTmr);


  mprintf("\nTiming Information --------------------------------------------------\n");
  mprintf("  I/O:          \t\t %7.3f\n", gk_getcputimer(IOTmr));
  mprintf("  Dual Creation:\t\t %7.3f\n", gk_getcputimer(DUALTmr));
  mprintf("**********************************************************************\n");

  }

  gk_free((void **)&elmnts, &xadj, &adjncy, &metype, &weights, &elms,  LTERM);
}


