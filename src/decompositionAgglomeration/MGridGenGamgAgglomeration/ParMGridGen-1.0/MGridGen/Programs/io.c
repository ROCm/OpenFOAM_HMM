/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * io.c
 *
 * This file contains routines related to I/O
 *
 * George Irene
 */

#include "mgridgen.h"


/*************************************************************************
* This function reads the spd matrix
**************************************************************************/
void ReadGraph(GraphType *graph, char *filename)
{
  int i, k, nvtxs, nedges;
  idxtype *xadj, *adjncy;
  realtype *vvol, *vsurf, *adjwgt;
  char line[MAXLINE+1], delim[] = " \t", *token;
  FILE *fpin;

  if ((fpin = fopen(filename, "r")) == NULL) {
    printf("Failed to open file %s\n", filename);
    exit(0);
  }

  do {
    fgets(line, MAXLINE, fpin);
  } while (line[0] == '%' && !feof(fpin));

  if (feof(fpin)) {
    graph->nvtxs = 0;
    return;
  }

  sscanf(line, "%d %d", &nvtxs, &nedges);
  nedges *= 2;

  graph->nvtxs = nvtxs;
  xadj   = graph->xadj   = idxsmalloc(nvtxs+1, 0, "ReadGraph: xadj");
  vvol   = graph->vvol   = realmalloc(nvtxs, "ReadGraph: vvol");
  vsurf  = graph->vsurf  = realsmalloc(nvtxs, 0.0, "ReadGraph: vsurf");
  adjncy = graph->adjncy = idxmalloc(nedges+1, "ReadGraph: adjncy");
  adjwgt = graph->adjwgt = realmalloc(nedges+1, "ReadGraph: adjwgt");

  /* Start reading the graph file */
  for (xadj[0]=0, k=0, i=0; i<nvtxs; i++) {
     do {
       fgets(line, MAXLINE, fpin);
     } while (line[0] == '%' && !feof(fpin));
     if (strlen(line) == MAXLINE) 
       errexit("\nBuffer for fgets not big enough!\n");

     /* Parse the string and get the arguments */
     token = strtok(line, delim);
     vvol[i] = atof(token);

     while ((token = strtok(NULL, delim))) {
       adjncy[k] = atoi(token) - 1;
       if (adjncy[k] == i)
         vsurf[i] = atof(strtok(NULL, delim));
       else
         adjwgt[k++] = atof(strtok(NULL, delim));
     }
     xadj[i+1] = k;
  }

  fclose(fpin);

  if (k != nedges)
    errexit("ReadGraph: Something wrong with the edges from input file %d %d",
            nedges, k);
}


/*************************************************************************
* Reads a matrix in mgridgen format and transforms it to a metis format
**************************************************************************/
void TransformGraph(char *filename)
{
  int i, length, nvtxs, nedges;
  idxtype adjncy;
  realtype vvol, vsurf, adjwgt;
  char line[MAXLINE+1], delim[] = " \t", *token;
  char *tfilename;
  FILE *fpin, *fpin2;

  length = strlen(filename);
  tfilename = (char *)malloc ((length +5) * sizeof(char));
  strcpy(tfilename, filename);
  strcat(tfilename, ".FORMETIS");

  if ((fpin = fopen(filename, "r")) == NULL) {
    printf("Failed to open file %s\n", filename);
    exit(0);
  }

  if ((fpin2 = fopen(tfilename, "w")) == NULL) {
    printf("Failed to open file %s\n", tfilename);
    exit(0);
  }

  do {
    fgets(line, MAXLINE, fpin);
  } while (line[0] == '%' && !feof(fpin));

  if (feof(fpin)) 
    return;

  sscanf(line, "%d %d", &nvtxs, &nedges);
  fprintf(fpin2, "%d %d\n", nvtxs, nedges);

  /* Start reading the graph file */
  for (i=0; i<nvtxs; i++) {
     do {
       fgets(line, MAXLINE, fpin);
     } while (line[0] == '%' && !feof(fpin));
     if (strlen(line) == MAXLINE) 
       errexit("\nBuffer for fgets not big enough!\n");

     /* Parse the string and get the arguments */
     token = strtok(line, delim);
     vvol = atof(token);

     while ((token = strtok(NULL, delim))) {
       adjncy = atoi(token);
       if (adjncy == i+1)
         vsurf = atof(strtok(NULL, delim));
       else {
         fprintf(fpin2, "%d ", adjncy);
         adjwgt = atof(strtok(NULL, delim));
       }
     }
     fprintf(fpin2, "\n");
  }

  fclose(fpin);
  fclose(fpin2);
}


/*************************************************************************
* This function writes out the partition vector
**************************************************************************/
void WritePartition(char *fname, idxtype *part, int n, int nparts)
{
  int i;
  char filename[256];
  FILE *fpout;

  sprintf(filename,"%s.part.%d",fname, nparts);

  if ((fpout = fopen(filename, "w")) == NULL) 
    errexit("Problems in opening the partition file: %s", filename);

  for (i=0; i<n; i++)
     fprintf(fpout,"%d\n",part[i]);

  fclose(fpout);
}


/*************************************************************************
* This function prints out the graph
**************************************************************************/
void PrintGraph(GraphType *graph)
{
  int i, j, nvtxs;
  idxtype *xadj, *adjncy;
  realtype *adjwgt;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  for (i=0; i<nvtxs; i++) {
     printf("%e %e ", graph->vvol[i], graph->vsurf[i]);
     for (j=xadj[i]; j<xadj[i+1]; j++)
        printf("%d %3f ", adjncy[j]+1, adjwgt[j]);
     printf("\n");
  }
}
