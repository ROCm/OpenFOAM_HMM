/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * parmgridgen.c
 * 
 * This file contains code for testing the adaptive partitioning routines
 *
 * George Irene
 */

#include "parmgridgen.h"


/***********************************************************************************
* This function is the testing routine for the adaptive multilevel partitioning code.
* It computes a partition from scratch, it then moves the graph and changes some
* of the vertex weights and then call the adaptive code.
************************************************************************************/
void TestParMGridGen(char *filename, int *options, int minsize, int maxsize, MPI_Comm comm)
{
  int i, nparts, npes, mype;
  MGridGraphType graph;
  idxtype *part;
  double tmr;


  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  MGridReadTestGraph(&graph, filename, comm);

  part = idxmalloc(graph.nvtxs, "TestParMGridGen: part");

  /*======================================================================
  / ParMETIS_AspectRatio
  /=======================================================================*/
  if (mype==0)
    printf("------------------------ PARAMETERS --------------------------------------\n");
  for (i=0; i<npes; i++)
    if (mype == i)
      printf("%s, Dim=%d [%2d %2d] CType=%d RType=%d Nvtxs=%d Nedges=%d\n", filename,
              options[OPTION_DIM], minsize, maxsize, options[OPTION_CTYPE],
              options[OPTION_RTYPE], graph.nvtxs, graph.nedges);

  cleartimer(tmr);
  MPI_Barrier(comm);
  starttimer(tmr);
  ParMGridGen(graph.vtxdist, graph.xadj, graph.vvol, graph.vsurf, graph.adjncy,
              graph.adjwgt, &nparts, minsize, maxsize, options, part, &comm);
  MPI_Barrier(comm);
  stoptimer(tmr);

  printf("Total Time = %lf\n", gettimer(tmr));

  WriteParallelPartition(filename, part, graph.vtxdist, nparts, mype, npes);

  IMfree(&graph.vtxdist, &graph.xadj, &graph.vvol, &graph.vsurf, &graph.vwgt,
         &graph.adjncy, &graph.adjwgt, &part, LTERM);
} 
