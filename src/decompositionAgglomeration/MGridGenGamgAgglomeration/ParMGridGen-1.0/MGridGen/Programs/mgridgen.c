/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * mgridgen.c
 *
 * This file contains the driving routine for MGRIDGEN
 *
 * George Irene
 */

#include "mgridgen.h"


/*************************************************************************
* Let the game begin
**************************************************************************/
int main(int argc, char *argv[])
{
  int options[10], nmoves, nparts, minsize, maxsize;
  idxtype *part;
  GraphType graph;
  char filename[256];
  double tmr;

  if (argc != 8) {
    printf("Usage: %s <GraphFile> <Dim> <CType> <RType> <minsize> <maxsize> <dbglvl>\n", argv[0]);
    printf("Where:\n");
    printf("\tDim:    \t2 -> 2-D Mesh\n");
    printf("\t        \t3 -> 3-D Mesh\n");
    printf("\tCType:  \t1 -> Random\n");
    printf("\t        \t2 -> HEM\n");
    printf("\t        \t3 -> Slow HEM\n");
    printf("\t        \t4 -> Slow Heaviest\n");
    printf("\tRType:  \t1 -> Aspect Ratio refinement\n");
    printf("\t        \t2 -> Weighted Aspect Ratio refinement\n");
    printf("\t        \t3 -> Surface cut refinement\n");
    printf("\t        \t4 -> Minimum Aspect Ratio & Average refinement\n");
    printf("\t        \t5 -> Minimum Aspect Ratio refinement\n");
    printf("\t        \t6 -> 4+2 Minimum & Weighted Aspect Ratio refinement\n");
    printf("\t        \t7 -> 5+2 Minimum & Weighted Aspect Ratio refinement\n");
    printf("\tminsize:\tA lower bound on the cell size (suggested)\n");
    printf("\tmaxsize:\tAn upper bound on the cell size (strict)\n");
    printf("--------------------------------------------------------------------\n");
    printf("Recomended usage: %s <GraphFile> <Dim> 4 6 ? ? 128\n", argv[0]);
    printf("--------------------------------------------------------------------\n");
    exit(0);
  }
    
  strcpy(filename, argv[1]);
  options[OPTION_DIM] = atoi(argv[2]);
  options[OPTION_CTYPE] = atoi(argv[3]);
  options[OPTION_RTYPE] = atoi(argv[4]);
  minsize = atoi(argv[5]);
  maxsize = atoi(argv[6]);
  options[OPTION_DBGLVL] = atoi(argv[7]);

  ReadGraph(&graph, filename);

  printf("------------------------ PARAMETERS --------------------------------------\n");
  printf("%s, Dim=%d [%2d %2d] CType=%d RType=%d Nvtxs=%d Nedges=%d\n", filename,
          options[OPTION_DIM], minsize, maxsize, options[OPTION_CTYPE], options[OPTION_RTYPE],
          graph.nvtxs, graph.xadj[graph.nvtxs]);

  part = idxmalloc(graph.nvtxs, "main: part");

  cleartimer(tmr);
  starttimer(tmr);
  MGridGen(graph.nvtxs, graph.xadj, graph.vvol, graph.vsurf, graph.adjncy, graph.adjwgt, 
           minsize, maxsize, options, &nmoves, &nparts, part);
  stoptimer(tmr);
  printf("Total Time: %lf\n", gettimer(tmr));

  WritePartition(filename, part, graph.nvtxs, nparts);

  IMfree(&graph.xadj, &graph.vvol, &graph.vsurf, &graph.adjncy, &graph.adjwgt, &part, LTERM);

  return(0);
}  
