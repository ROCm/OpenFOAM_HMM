/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * metis.c
 *
 * This file contains the driving routine for multilevel method
 *
 * Started 8/28/94
 * George
 *
 * $Id: metis.c,v 1.2 2002/08/10 06:02:53 karypis Exp $
 *
 */

#include <metisbin.h>



/*************************************************************************
* Let the game begin
**************************************************************************/
int main(int argc, char *argv[])
{
  idxtype i, nparts, OpType, options[10], nbytes;
  idxtype *part, *perm, *iperm, *sizes;
  GraphType graph;
  char filename[256];
  idxtype numflag = 0, wgtflag = 0, edgecut;


  if (argc != 11) {
    mprintf("Usage: %s <GraphFile> <Nparts> <Mtype> <Rtype> <IPtype> <Oflags> <Pfactor> <Nseps> <OPtype> <Options> \n",argv[0]);
    exit(0);
  }
    
  strcpy(filename, argv[1]);
  nparts = atoi(argv[2]);
  options[OPTION_CTYPE] = atoi(argv[3]);
  options[OPTION_RTYPE] = atoi(argv[4]);
  options[OPTION_ITYPE] = atoi(argv[5]);
  options[OPTION_OFLAGS] = atoi(argv[6]);
  options[OPTION_PFACTOR] = atoi(argv[7]);
  options[OPTION_NSEPS] = atoi(argv[8]);
  OpType = atoi(argv[9]); 
  options[OPTION_DBGLVL] = atoi(argv[10]);


  ReadGraph(&graph, filename, &wgtflag);
  if (graph.nvtxs <= 0) {
    mprintf("Empty graph. Nothing to do.\n");
    exit(0);
  }
  mprintf("Partitioning a graph with %D vertices and %D edges\n", graph.nvtxs, graph.nedges/2);

  METIS_EstimateMemory(&graph.nvtxs, graph.xadj, graph.adjncy, &numflag, &OpType, &nbytes);
  mprintf("Metis will need %D Mbytes\n", nbytes/(1024*1024));

  part = perm = iperm = NULL;

  options[0] = 1;
  switch (OpType) {
    case OP_PMETIS:
      mprintf("Recursive Partitioning... ------------------------------------------\n");
      part = idxmalloc(graph.nvtxs, "main: part");

      METIS_PartGraphRecursive(&graph.nvtxs, graph.xadj, graph.adjncy, graph.vwgt, graph.adjwgt, 
                               &wgtflag, &numflag, &nparts, options, &edgecut, part);

      IFSET(options[OPTION_DBGLVL], DBG_OUTPUT, WritePartition(filename, part, graph.nvtxs, nparts)); 

      mprintf("  %D-way Edge-Cut: %7D\n", nparts, edgecut);
      ComputePartitionInfo(&graph, nparts, part);

      gk_free((void **)&part, LTERM);
      break;
    case OP_KMETIS:
      mprintf("K-way Partitioning... -----------------------------------------------\n");
      part = idxmalloc(graph.nvtxs, "main: part");

      METIS_PartGraphKway(&graph.nvtxs, graph.xadj, graph.adjncy, graph.vwgt, graph.adjwgt, 
                          &wgtflag, &numflag, &nparts, options, &edgecut, part);

      IFSET(options[OPTION_DBGLVL], DBG_OUTPUT, WritePartition(filename, part, graph.nvtxs, nparts)); 

      mprintf("  %D-way Edge-Cut: %7D\n", nparts, edgecut);
      ComputePartitionInfo(&graph, nparts, part);

      gk_free((void **)&part, LTERM);
      break;
    case OP_OEMETIS:
      gk_free((void **)&graph.vwgt, &graph.adjwgt, LTERM);

      mprintf("Edge-based Nested Dissection Ordering... ----------------------------\n");
      perm = idxmalloc(graph.nvtxs, "main: perm");
      iperm = idxmalloc(graph.nvtxs, "main: iperm");

      METIS_EdgeND(&graph.nvtxs, graph.xadj, graph.adjncy, &numflag, options, perm, iperm);

      IFSET(options[OPTION_DBGLVL], DBG_OUTPUT, WritePartition(filename, iperm, graph.nvtxs, 0)); 

      ComputeFillIn(&graph, iperm);

      gk_free((void **)&perm, &iperm, LTERM);
      break;
    case OP_ONMETIS:
      gk_free((void **)&graph.vwgt, &graph.adjwgt, LTERM);

      mprintf("Node-based Nested Dissection Ordering... ----------------------------\n");
      perm = idxmalloc(graph.nvtxs, "main: perm");
      iperm = idxmalloc(graph.nvtxs, "main: iperm");

      METIS_NodeND(&graph.nvtxs, graph.xadj, graph.adjncy, &numflag, options, perm, iperm);

      IFSET(options[OPTION_DBGLVL], DBG_OUTPUT, WritePartition(filename, iperm, graph.nvtxs, 0)); 

      ComputeFillIn(&graph, iperm);

      gk_free((void **)&perm, &iperm, LTERM);
      break;
    case OP_ONWMETIS:
      gk_free((void **)&graph.adjwgt, LTERM);

      mprintf("WNode-based Nested Dissection Ordering... ---------------------------\n");
      perm = idxmalloc(graph.nvtxs, "main: perm");
      iperm = idxmalloc(graph.nvtxs, "main: iperm");

      METIS_NodeWND(&graph.nvtxs, graph.xadj, graph.adjncy, graph.vwgt, &numflag, options, perm, iperm);

      IFSET(options[OPTION_DBGLVL], DBG_OUTPUT, WritePartition(filename, iperm, graph.nvtxs, 0)); 

      ComputeFillIn(&graph, iperm);

      gk_free((void **)&perm, &iperm, LTERM);
      break;
    case 6:
      gk_free((void **)&graph.vwgt, &graph.adjwgt, LTERM);

      mprintf("Node-based Nested Dissection Ordering... ----------------------------\n");
      perm = idxmalloc(graph.nvtxs, "main: perm");
      iperm = idxmalloc(graph.nvtxs, "main: iperm");
      sizes = idxmalloc(2*nparts, "main: sizes");

      METIS_NodeNDP(graph.nvtxs, graph.xadj, graph.adjncy, nparts, options, perm, iperm, sizes);

      IFSET(options[OPTION_DBGLVL], DBG_OUTPUT, WritePartition(filename, iperm, graph.nvtxs, 0)); 

      ComputeFillIn(&graph, iperm);

      for (i=0; i<2*nparts-1; i++)
        mprintf("%D ", sizes[i]);
      mprintf("\n");

      gk_free((void **)&perm, &iperm, &sizes, LTERM);
      break;
    case 7:
      mprintf("K-way Vol Partitioning... -------------------------------------------\n");
      part = idxmalloc(graph.nvtxs, "main: part");

      METIS_PartGraphVKway(&graph.nvtxs, graph.xadj, graph.adjncy, graph.vwgt, NULL, 
            &wgtflag, &numflag, &nparts, options, &edgecut, part);

      IFSET(options[OPTION_DBGLVL], DBG_OUTPUT, WritePartition(filename, part, graph.nvtxs, nparts)); 

      mprintf("  %D-way Volume: %7D\n", nparts, edgecut);
      ComputePartitionInfo(&graph, nparts, part);

      gk_free((void **)&part, LTERM);
      break;
    case 9:
      mprintf("K-way Partitioning (with vwgts)... ----------------------------------\n");
      part = idxmalloc(graph.nvtxs, "main: part");
      graph.vwgt = idxmalloc(graph.nvtxs, "main: graph.vwgt");
      for (i=0; i<graph.nvtxs; i++)
        graph.vwgt[i] = graph.xadj[i+1]-graph.xadj[i]+1;
      wgtflag = 2;

      METIS_PartGraphKway(&graph.nvtxs, graph.xadj, graph.adjncy, graph.vwgt, graph.adjwgt, 
                          &wgtflag, &numflag, &nparts, options, &edgecut, part);

      IFSET(options[OPTION_DBGLVL], DBG_OUTPUT, WritePartition(filename, part, graph.nvtxs, nparts)); 

      mprintf("  %D-way Edge-Cut: %7D\n", nparts, edgecut);
      ComputePartitionInfo(&graph, nparts, part);

      gk_free((void **)&part, LTERM);
      break;
    case 10:
      break;
    default:
      errexit("Unknown");
  }

  gk_free((void **)&graph.xadj, &graph.adjncy, &graph.vwgt, &graph.adjwgt, LTERM);
}  


