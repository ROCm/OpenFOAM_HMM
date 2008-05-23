/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * debug.c
 *
 * This file contains various functions that are used to display debuging 
 * information
 *
 * George Irene
 */

#include "parmetis.h"


/*************************************************************************
* This function prints a vector stored in each processor 
**************************************************************************/
void PrintVector2(CtrlType *ctrl, int n, int first, idxtype *vec, char *title)
{
  int i, penum;

  for (penum=0; penum<ctrl->npes; penum++) {
    if (ctrl->mype == penum) {
      if (ctrl->mype == 0)
        printf("%s\n", title);
      printf("\t%3d. ", ctrl->mype);
      for (i=0; i<n; i++)
        printf("[%d %d.%hd] ", first+i, (vec[i]>=KEEP_BIT ? 1 : 0), (vec[i]>=KEEP_BIT ? vec[i]-KEEP_BIT : vec[i]));
      printf("\n");
      fflush(stdout);
    }
    MPI_Barrier(ctrl->comm);
  }
}
