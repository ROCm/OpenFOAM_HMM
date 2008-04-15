/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * main.c
 * 
 * This file contains code for testing the adaptive partitioning routines
 *
 * George Irene
 */

#include "parmgridgen.h"


/*************************************************************************
* Let the game begin
**************************************************************************/
int main(int argc, char *argv[])
{
  int mype, npes;
  MPI_Comm comm;
  int minsize, maxsize;
  int options[10];
  char filename[15];

  MPI_Init(&argc, &argv);
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  sprintf(filename,"OUTPUT%d.%d",npes,mype);
  freopen(filename, "a", stdout);

  if (argc != 8) {
    if (mype == 0) {
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
    }

    MPI_Finalize();
    exit(0);
  }

  options[OPTION_DIM] = atoi(argv[2]);
  options[OPTION_CTYPE] = atoi(argv[3]);
  options[OPTION_RTYPE] = atoi(argv[4]);
  minsize = atoi(argv[5]);
  maxsize = atoi(argv[6]);
  options[OPTION_DBGLVL] = atoi(argv[7]);
  TestParMGridGen(argv[1], options, minsize, maxsize, comm); 

  MPI_Comm_free(&comm);

  MPI_Finalize();
  return(0);
}
