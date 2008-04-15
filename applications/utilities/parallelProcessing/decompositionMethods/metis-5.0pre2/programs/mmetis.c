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
 * $Id: mmetis.c,v 1.2 2002/08/10 06:02:53 karypis Exp $
 *
 */

#include <metisbin.h>



/*************************************************************************
* Let the game begin
**************************************************************************/
int main(int argc, char *argv[])
{
  idxtype i, j, nparts, OpType, options[10], nbytes;
  idxtype *part, *perm, *iperm, *sizes;
  GraphType graph;
  char filename[256];
  idxtype numflag = 0, wgtflag = 0, edgecut;
  idxtype *mvwgt, *vwgt, *vsize;
  float rubvec[MAXNCON];


  if (argc < 8) {
    mprintf("Usage: %s <GraphFile> <Nparts> <Mtype> <Rtype> <IPtype> <OpType> <Options> \n",argv[0]);
    exit(0);
  }
    
  strcpy(filename, argv[1]);
  nparts = atoi(argv[2]);
  options[OPTION_CTYPE] = atoi(argv[3]);
  options[OPTION_RTYPE] = atoi(argv[4]);
  options[OPTION_ITYPE] = atoi(argv[5]);
  OpType = atoi(argv[6]); 
  options[OPTION_DBGLVL] = atoi(argv[7]);


  ReadGraph(&graph, filename, &wgtflag);
  if (graph.nvtxs <= 0) {
    mprintf("Empty graph. Nothing to do.\n");
    exit(0);
  }
  mprintf("Partitioning a graph with %D vertices and %D edges. Constraints: %D\n", graph.nvtxs, graph.nedges/2, graph.ncon);

  part = perm = iperm = NULL;
  vsize = NULL;

  options[0] = 1;
  switch (OpType) {
    case OP_PMETIS:
      mprintf("Recursive Partitioning... ------------------------------------------\n");
      part = idxmalloc(graph.nvtxs, "main: part");

      METIS_mCPartGraphRecursive(&graph.nvtxs, &graph.ncon, graph.xadj, graph.adjncy, 
         graph.vwgt, graph.adjwgt, &wgtflag, &numflag, &nparts, options, &edgecut, part);

      IFSET(options[OPTION_DBGLVL], DBG_OUTPUT, WritePartition(filename, part, graph.nvtxs, nparts)); 

      mprintf("  %D-way Edge-Cut: %7D\n", nparts, edgecut);
      ComputePartitionInfo(&graph, nparts, part); 

      gk_free((void **)&part, LTERM);
      break;
    case OP_KMETIS:
      mprintf("K-way Partitioning... ----------------------------------------------\n");
      part = idxmalloc(graph.nvtxs, "main: part");

      if (argc != 8+graph.ncon) 
        errexit("You must supply %d ub constraints!\n", graph.ncon);

      for (i=0; i<graph.ncon; i++)
        rubvec[i] = atof(argv[8+i]);

      METIS_mCPartGraphKway(&graph.nvtxs, &graph.ncon, graph.xadj, graph.adjncy, graph.vwgt, 
         graph.adjwgt, &wgtflag, &numflag, &nparts, rubvec, options, &edgecut, part);

      IFSET(options[OPTION_DBGLVL], DBG_OUTPUT, WritePartition(filename, part, graph.nvtxs, nparts)); 

      mprintf("  %D-way Edge-Cut: %7D \tUBVec: ", nparts, edgecut);
      for (i=0; i<graph.ncon; i++)
        mprintf("%.3f ", rubvec[i]);
      mprintf("\n");

      ComputePartitionInfo(&graph, nparts, part); 

      if (options[OPTION_DBGLVL]&1024) {
        /* Partition each objective separately and see the results */
        vwgt = idxmalloc(graph.nvtxs, "vwgt");

        for (j=0; j<graph.ncon; j++) {
          for (i=0; i<graph.nvtxs; i++)
            vwgt[i] = graph.vwgt[i*graph.ncon+j];

          options[0] = 0;
          METIS_PartGraphKway(&graph.nvtxs, graph.xadj, graph.adjncy, vwgt, graph.adjwgt, 
               &wgtflag, &numflag, &nparts, options, &edgecut, part);

          mprintf("Partitioning using constrain %D ------------------------------------\n", j);
          ComputePartitionInfo(&graph, nparts, part);
        }
        gk_free((void **)&vwgt, LTERM);
      }

      gk_free((void **)&part, LTERM);
      break;
    case 3:
      mprintf("Recursive Partitioning... -----------------------------------------\n");
      part = idxmalloc(graph.nvtxs, "main: part");

      if (argc != 8+graph.ncon) 
        errexit("You must supply %d ub constraints!\n", graph.ncon);

      for (i=0; i<graph.ncon; i++)
        rubvec[i] = atof(argv[8+i]);

      METIS_mCHPartGraphRecursive(&graph.nvtxs, &graph.ncon, graph.xadj, graph.adjncy, 
         graph.vwgt, graph.adjwgt, &wgtflag, &numflag, &nparts, rubvec, options, &edgecut, part);

      IFSET(options[OPTION_DBGLVL], DBG_OUTPUT, WritePartition(filename, part, graph.nvtxs, nparts)); 

      mprintf("  %D-way Edge-Cut: %7D \tUBVec: ", nparts, edgecut);
      for (i=0; i<graph.ncon; i++)
        mprintf("%.3f ", rubvec[i]);
      mprintf("\n");

      ComputePartitionInfo(&graph, nparts, part); 

      gk_free((void **)&part, LTERM);
      break;
    default:
      errexit("Unknown");
  }

  gk_free((void **)&graph.xadj, &graph.adjncy, &graph.vwgt, &graph.adjwgt, LTERM);
}  


