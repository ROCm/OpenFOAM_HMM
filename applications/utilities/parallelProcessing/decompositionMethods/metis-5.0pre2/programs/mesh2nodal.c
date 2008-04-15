/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mesh2nodal.c
 *
 * This file reads in the element node connectivity array of a mesh and writes
 * out its dual in the format suitable for Metis.
 *
 * Started 9/29/97
 * George
 *
 * $Id: mesh2nodal.c,v 1.2 2002/08/10 06:02:53 karypis Exp $
 *
 */

#include <metisbin.h>



/*************************************************************************
* Let the game begin
**************************************************************************/
int main(int argc, char *argv[])
{
  idxtype i, j, ne, nn, etype, mtype, numflag=0;
  idxtype *elmnts, *xadj, *adjncy, *metype, *weights;
  double IOTmr, DUALTmr;
  char fileout[256], etypestr[5][5] = {"TRI", "TET", "HEX", "QUAD", "LINE"};

  if (argc != 2) {
    mprintf("Usage: %s <meshfile>\n",argv[0]);
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
  mprintf("Forming Nodal Graph... ----------------------------------------------\n");

  xadj = idxmalloc(nn+1, "main: xadj");
  adjncy = idxmalloc(20*nn, "main: adjncy");

  gk_startcputimer(DUALTmr);
  METIS_MeshToNodal(&ne, &nn, elmnts, &etype, &numflag, xadj, adjncy);
  gk_stopcputimer(DUALTmr);

  mprintf("  Nodal Information: #Vertices: %D, #Edges: %D\n", nn, xadj[nn]/2);

  msprintf(fileout, "%s.ngraph", argv[1]);
  gk_startcputimer(IOTmr);
  WriteGraph(fileout, nn, xadj, adjncy);
  gk_stopcputimer(IOTmr);


  mprintf("\nTiming Information --------------------------------------------------\n");
  mprintf("  I/O:          \t\t %7.3f\n", gk_getcputimer(IOTmr));
  mprintf("  Nodal Creation:\t\t %7.3f\n", gk_getcputimer(DUALTmr));
  mprintf("**********************************************************************\n");


  }

  else{

  gk_clearcputimer(IOTmr);
  gk_clearcputimer(DUALTmr);

  gk_startcputimer(IOTmr);
  if(mtype==0)
     elmnts = ReadMixedMesh(argv[1], &ne, &nn, metype);
  else
     elmnts = ReadMixedMeshWgt(argv[1], &ne, &nn, metype, weights);
  gk_stopcputimer(IOTmr);


  mprintf("**********************************************************************\n");
  mprintf("%s", METISTITLE);
  mprintf("Mesh Information ----------------------------------------------------\n");
  mprintf("  Name: %s, #Elements: %D, #Nodes: %D, Etype: %s\n\n", argv[1], ne, nn, "Mixed");
  mprintf("Forming Nodal Graph... ----------------------------------------------\n");

  xadj = idxmalloc(nn+1, "main: xadj");
  adjncy = idxmalloc(20*nn, "main: adjncy");

  gk_startcputimer(DUALTmr);
  METIS_MixedMeshToNodal(&ne, &nn, elmnts, metype, &numflag, xadj, adjncy);
  gk_stopcputimer(DUALTmr);

  mprintf("  Nodal Information: #Vertices: %D, #Edges: %D\n", nn, xadj[nn]/2);

  msprintf(fileout, "%s.ngraph", argv[1]);
  gk_startcputimer(IOTmr);
  WriteGraph(fileout, nn, xadj, adjncy);
  gk_stopcputimer(IOTmr);


  mprintf("\nTiming Information --------------------------------------------------\n");
  mprintf("  I/O:          \t\t %7.3f\n", gk_getcputimer(IOTmr));
  mprintf("  Nodal Creation:\t\t %7.3f\n", gk_getcputimer(DUALTmr));
  mprintf("**********************************************************************\n");


  }
  
  gk_free((void **)&elmnts, &xadj, &adjncy, &metype, &weights, LTERM);

}


