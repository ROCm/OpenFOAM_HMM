/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * partdmesh.c
 *
 * This file reads in the element node connectivity array of a mesh and 
 * partitions both the elements and the nodes using KMETIS on the dual graph.
 *
 * Started 9/29/97
 * George
 *
 * $Id: partdmesh.c,v 1.2 2002/08/10 06:02:54 karypis Exp $
 *
 */

#include <metisbin.h>



/*************************************************************************
* Let the game begin
**************************************************************************/
int main(int argc, char *argv[])
{
  idxtype i, j, ne, nn, etype, mtype, numflag=0, nparts, edgecut, custom=0;
  idxtype *elmnts, *epart, *npart, *metype, *conmat, *weights;
  double IOTmr, DUALTmr;
  char etypestr[5][5] = {"TRI", "TET", "HEX", "QUAD", "LINE"};
  GraphType graph;

  if (argc < 3) {
    mprintf("Usage: %s <meshfile> <nparts> [confile]\n",argv[0]);
    exit(0);
  }

  nparts = atoi(argv[2]);
  if (nparts < 2) {
    mprintf("nparts must be greater than one.\n");
    exit(0);
  }
  mtype=MeshType(argv[1]);
  ne=MixedElements(argv[1]);
  metype = idxmalloc(ne, "main: metype");
  weights = idxmalloc(ne, "main: weights");
 
  gk_clearcputimer(IOTmr);
  gk_clearcputimer(DUALTmr);

  gk_startcputimer(IOTmr);

  if(mtype==1)
       elmnts = ReadMesh(argv[1], &ne, &nn, &etype);
  else if(mtype==3)
       elmnts = ReadMeshWgt(argv[1], &ne, &nn, &etype, weights);
  else if(mtype==0)
       elmnts = ReadMixedMesh(argv[1], &ne, &nn, metype);
  else
       elmnts = ReadMixedMeshWgt(argv[1], &ne, &nn, metype, weights);


  if (argc==4){
  conmat = ReadMgcnums(argv[2]);
  custom=1;
  }

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
  mprintf("Partitioning Dual Graph... ------------------------------------------\n");


  gk_startcputimer(DUALTmr);
  if (mtype==1)
  METIS_PartMeshDual(&ne, &nn, elmnts, &etype, &numflag, &nparts, &edgecut, epart, npart, 0, NULL);
  else if (mtype==3)
  METIS_PartMeshDual(&ne, &nn, elmnts, &etype, &numflag, &nparts, &edgecut, epart, npart, 2, weights);
  else if (mtype==0)
  METIS_PartMixedMeshDual(&ne, &nn, elmnts, metype, &numflag, &nparts, &edgecut, epart, npart, conmat, custom, 0, NULL);
  else 
  METIS_PartMixedMeshDual(&ne, &nn, elmnts, metype, &numflag, &nparts, &edgecut, epart, npart, conmat, custom, 2, weights);
  
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
  graph.nvtxs = nn;
  graph.xadj = idxmalloc(nn+1, "xadj");
  graph.vwgt = idxsmalloc(nn, 1, "vwgt");
  graph.adjncy = idxmalloc(20*nn, "adjncy");
  graph.adjwgt = idxsmalloc(20*nn, 1, "adjncy");

  METIS_MeshToNodal(&ne, &nn, elmnts, &etype, &numflag, graph.xadj, graph.adjncy);

  ComputePartitionInfo(&graph, nparts, npart);

  gk_free((void **)&graph.xadj, &graph.adjncy, &graph.vwgt, &graph.adjwgt, LTERM);
*/

  gk_free((void **)&elmnts, &epart, &npart, &metype, &weights, LTERM);

}


