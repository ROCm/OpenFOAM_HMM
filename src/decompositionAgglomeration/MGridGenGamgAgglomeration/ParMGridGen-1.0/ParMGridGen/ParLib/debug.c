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

#include "parmgridgen.h"


/*************************************************************************
* This function prints the local portion of the graph stored at each 
* processor
**************************************************************************/
void PrintGraph(MGridCtrlType *ctrl, MGridGraphType *graph)
{
  int i, j, penum;
  int firstvtx;

  MPI_Barrier(ctrl->comm);

  firstvtx = graph->vtxdist[ctrl->mype];

  for (penum=0; penum<ctrl->npes; penum++) {
    if (ctrl->mype == penum) {
      printf("\t%d", penum);
      for (i=0; i<graph->nvtxs; i++) {
        if (i==0)
          printf("\t%2d %2d\t", firstvtx+i, graph->vwgt[i]);
        else
          printf("\t\t%2d %2d\t", firstvtx+i, graph->vwgt[i]);
        for (j=graph->xadj[i]; j<graph->xadj[i+1]; j++)
          printf("[%d %f] ", graph->adjncy[j], graph->adjwgt[j]);
        printf("\n");
      }
      fflush(stdout);
    }
    MPI_Barrier(ctrl->comm);
  }
}


/*************************************************************************
* This function prints the local portion of the graph stored at each 
* processor along with degree information during refinement
**************************************************************************/
void PrintGraph2(MGridCtrlType *ctrl, MGridGraphType *graph)
{
  int i, j, penum;
  int firstvtx;

  MPI_Barrier(ctrl->comm);

  firstvtx = graph->vtxdist[ctrl->mype];

  for (penum=0; penum<ctrl->npes; penum++) {
    if (ctrl->mype == penum) {
      printf("\t%d", penum);
      for (i=0; i<graph->nvtxs; i++) {
        if (i==0)
          printf("\t%2d %2d [%d %f %f]\t", firstvtx+i, graph->vwgt[i], graph->where[i], graph->rinfo[i].id, graph->rinfo[i].ed);
        else
          printf("\t\t%2d %2d [%d %f %f]\t", firstvtx+i, graph->vwgt[i], graph->where[i], graph->rinfo[i].id, graph->rinfo[i].ed);
        for (j=graph->xadj[i]; j<graph->xadj[i+1]; j++)
          printf("[%d %f] ", graph->adjncy[j], graph->adjwgt[j]);
        printf("\n");
      }
      fflush(stdout);
    }
    MPI_Barrier(ctrl->comm);
  }
}

/*************************************************************************
* This function writes a graph in the format used by serial METIS
**************************************************************************/
void WriteMetisGraph(int nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, realtype *adjwgt)
{
  int i, j;
  FILE *fp;

  fp = fopen("test.graph", "w");

  fprintf(fp, "%d %d 11", nvtxs, xadj[nvtxs]/2);
  for (i=0; i<nvtxs; i++) {
    fprintf(fp, "\n%d ", vwgt[i]);
    for (j=xadj[i]; j<xadj[i+1]; j++)
      fprintf(fp, " %d %f", adjncy[j]+1, adjwgt[j]);
  }
  fclose(fp);
}

void WriteMetisGraph2(int nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, realtype *adjwgt, MPI_Comm comm)
{
   int i, j;
   int mype, npes;
   char filename[20];
   FILE *fp;

   MPI_Comm_size(comm, &npes);
   MPI_Comm_rank(comm, &mype);

   sprintf(filename,"test.graph%d.%d",npes,mype);
   fp = fopen(filename, "w");

   fprintf(fp, "%d %d", nvtxs, xadj[nvtxs]/2);

   for (i=0; i<nvtxs; i++) {
      fprintf(fp, "\n%d ", vwgt[i]);
      for (j=xadj[i]; j<xadj[i+1]; j++)
      fprintf(fp, " %d %f", adjncy[j]+1, adjwgt[j]);
  }
  fclose(fp);
}
