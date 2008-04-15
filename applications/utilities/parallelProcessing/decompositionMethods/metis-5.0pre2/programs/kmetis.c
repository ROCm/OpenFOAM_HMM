/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * kmetis.c
 *
 * This file contains the driving routine for kmetis
 *
 * Started 8/28/94
 * George
 *
 * $Id: kmetis.c,v 1.6 2003/07/31 16:15:50 karypis Exp $
 *
 */

#include <metisbin.h>



/*************************************************************************
* Let the game begin
**************************************************************************/
int main(int argc, char *argv[])
{
  idxtype i, nparts, options[10];
  idxtype *part;
  float rubvec[MAXNCON], lbvec[MAXNCON];
  GraphType graph;
  char filename[256];
  idxtype numflag = 0, wgtflag = 0, edgecut;
  double TOTALTmr, METISTmr, IOTmr;

  if (argc != 3) {
    mprintf("Usage: %s <GraphFile> <Nparts>\n",argv[0]);
    exit(0);
  }
    
  strcpy(filename, argv[1]);
  nparts = atoi(argv[2]);

  if (nparts < 2) {
    mprintf("The number of partitions should be greater than 1!\n");
    exit(0);
  }

  gk_clearcputimer(TOTALTmr);
  gk_clearcputimer(METISTmr);
  gk_clearcputimer(IOTmr);

  gk_startcputimer(TOTALTmr);
  gk_startcputimer(IOTmr);
  ReadGraph(&graph, filename, &wgtflag);

  /* The following is for debuging empty graphs... 
  graph.nedges = 0;
  idxset(graph.nvtxs+1, 0, graph.xadj);
  */

  if (graph.nvtxs <= 0) {
    mprintf("Empty graph. Nothing to do.\n");
    exit(0);
  }
  gk_stopcputimer(IOTmr);

  mprintf("**********************************************************************\n");
  mprintf("%s", METISTITLE);
  mprintf("Graph Information ---------------------------------------------------\n");
  mprintf("  Name: %s, #Vertices: %D, #Edges: %D, #Parts: %D\n", filename, graph.nvtxs, graph.nedges/2, nparts);
  if (graph.ncon > 1)
    mprintf("  Balancing Constraints: %D\n", graph.ncon);
  mprintf("\nK-way Partitioning... -----------------------------------------------\n");

  part = idxmalloc(graph.nvtxs, "main: part");
  options[0] = 0;

  gk_startcputimer(METISTmr);
  if (graph.ncon == 1) {
    METIS_PartGraphKway(&graph.nvtxs, graph.xadj, graph.adjncy, graph.vwgt, graph.adjwgt, 
          &wgtflag, &numflag, &nparts, options, &edgecut, part);
  }
  else {
    for (i=0; i<graph.ncon; i++)
      rubvec[i] = HORIZONTAL_IMBALANCE;

    METIS_mCPartGraphKway(&graph.nvtxs, &graph.ncon, graph.xadj, graph.adjncy, graph.vwgt, 
          graph.adjwgt, &wgtflag, &numflag, &nparts, rubvec, options, &edgecut, part);
  }
  gk_stopcputimer(METISTmr);

  ComputePartitionBalance(&graph, nparts, part, lbvec);

  mprintf("  %D-way Edge-Cut: %7D, Balance: ", nparts, edgecut);
  for (i=0; i<graph.ncon; i++)
    mprintf("%5.2f ", lbvec[i]);
  mprintf("\n");

  gk_startcputimer(IOTmr);
  WritePartition(filename, part, graph.nvtxs, nparts); 
  gk_stopcputimer(IOTmr);
  gk_stopcputimer(TOTALTmr);

  mprintf("\nTiming Information --------------------------------------------------\n");
  mprintf("  I/O:          \t\t %7.3f\n", gk_getcputimer(IOTmr));
  mprintf("  Partitioning: \t\t %7.3f   (KMETIS time)\n", gk_getcputimer(METISTmr));
  mprintf("  Total:        \t\t %7.3f\n", gk_getcputimer(TOTALTmr));
  mprintf("**********************************************************************\n");


  gk_free((void **)&graph.xadj, &graph.adjncy, &graph.vwgt, &graph.adjwgt, &part, LTERM);
}  


