/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * partnmesh.c
 *
 * This file reads in the element node connectivity array of a mesh and 
 * partitions both the elements and the nodes using KMETIS on the dual graph.
 *
 * Started 9/29/97
 * George
 *
 * $Id: partnmesh.c,v 1.2 2002/08/10 06:02:54 karypis Exp $
 *
 */

#include <metisbin.h>



/*************************************************************************
* Let the game begin
**************************************************************************/
int main(int argc, char *argv[])
{
  idxtype i, j, ne, nn, etype, mtype, numflag=0, nparts, edgecut;
  idxtype *elmnts, *epart, *npart, *metype, *weights;
  double IOTmr, DUALTmr;
  char etypestr[5][5] = {"TRI", "TET", "HEX", "QUAD", "LINE"};
  GraphType graph;

  if (argc != 3) {
    mprintf("Usage: %s <meshfile> <nparts>\n",argv[0]);
    exit(0);
  }



  nparts = atoi(argv[2]);
  if (nparts < 2) {
    mprintf("nparts must be greater than one.\n");
    exit(0);
  }
   
  gk_clearcputimer(IOTmr);
  gk_clearcputimer(DUALTmr);
 
  mtype=MeshType(argv[1]);
  ne=MixedElements(argv[1]);
  metype = idxmalloc(ne, "main: metype");
  weights = idxmalloc(ne, "main: weights");

  gk_startcputimer(IOTmr);
 
  if(mtype==1)
       elmnts = ReadMesh(argv[1], &ne, &nn, &etype);
  else if(mtype==3)
       elmnts = ReadMeshWgt(argv[1], &ne, &nn, &etype, weights);
  else if(mtype==0)
       elmnts = ReadMixedMesh(argv[1], &ne, &nn, metype);
  else
       elmnts = ReadMixedMeshWgt(argv[1], &ne, &nn, metype, weights);

  gk_stopcputimer(IOTmr);

  epart = idxmalloc(ne, "main: epart");
  npart = idxmalloc(nn, "main: npart");

  mprintf("**********************************************************************\n");
  mprintf("%s", METISTITLE);
  mprintf("Mesh Information ----------------------------------------------------\n");
  if (mtype==1)
  mprintf("  Name: %s, #Elements: %D, #Nodes: %D, Etype: %s\n\n", argv[1], ne, nn, etypestr[etype-1]);
  else
  mprintf("  Name: %s, #Elements: %D, #Nodes: %D, Etype: %s\n\n", argv[1], ne, nn, "Mixed");
  mprintf("Partitioning Nodal Graph... -----------------------------------------\n");


  gk_startcputimer(DUALTmr);
  
  if (mtype==1 || mtype==3)
  METIS_PartMeshNodal(&ne, &nn, elmnts, &etype, &numflag, &nparts, &edgecut, epart, npart);
  else 
  METIS_PartMixedMeshNodal(&ne, &nn, elmnts, metype, &numflag, &nparts, &edgecut, epart, npart);
  
  gk_stopcputimer(DUALTmr);

  mprintf("  %D-way Edge-Cut: %7D, Balance: %5.2f\n", nparts, edgecut, ComputeElementBalance(ne, nparts, epart));

  gk_startcputimer(IOTmr);
  WriteMeshPartition(argv[1], nparts, ne, epart, nn, npart);
  gk_stopcputimer(IOTmr);


  mprintf("\nTiming Information --------------------------------------------------\n");
  mprintf("  I/O:          \t\t %7.3f\n", gk_getcputimer(IOTmr));
  mprintf("  Partitioning: \t\t %7.3f\n", gk_getcputimer(DUALTmr));
  mprintf("**********************************************************************\n");

/*
  graph.nvtxs = ne;
  graph.xadj = idxmalloc(ne+1, "xadj");
  graph.vwgt = idxsmalloc(ne, 1, "vwgt");
  graph.adjncy = idxmalloc(10*ne, "adjncy");
  graph.adjwgt = idxsmalloc(10*ne, 1, "adjncy");

  METIS_MeshToDual(&ne, &nn, elmnts, &etype, &numflag, graph.xadj, graph.adjncy);

  ComputePartitionInfo(&graph, nparts, epart);

  gk_free((void **)&graph.xadj, &graph.adjncy, &graph.vwgt, &graph.adjwgt, LTERM);
*/

  gk_free((void **)&elmnts, &epart, &npart, &metype, &weights,  LTERM);

}


