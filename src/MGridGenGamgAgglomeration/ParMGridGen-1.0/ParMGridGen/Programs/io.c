/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * tio.c
 *
 * This file contains routines related to I/O
 *
 * George Irene
 */

#include "parmgridgen.h"


/*************************************************************************
* This function reads the CSR matrix
**************************************************************************/
void MGridReadTestGraph(MGridGraphType *graph, char *filename, MPI_Comm comm)
{
  int i, k, l, npes, mype;
  int nvtxs, penum, snvtxs;
  idxtype *gxadj, *gadjncy;  
  idxtype *vtxdist, *sxadj, *ssize;
  realtype *gadjwgt, *gvsurf, *gvvol;
  MPI_Status status;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  vtxdist = graph->vtxdist = idxsmalloc(npes+1, 0, "MGridReadTestGraph: vtxdist");

  if (mype == 0) {
    ssize = idxsmalloc(npes, 0, "MGridReadTestGraph: ssize");

    ReadMGridGraph(filename, &nvtxs, &gxadj, &gadjncy, &gadjwgt, &gvvol, &gvsurf);

    printf("Nvtxs: %d, Nedges: %d\n", nvtxs, gxadj[nvtxs]);

    /* Construct vtxdist and send it to all the processors */
    vtxdist[0] = 0;
    for (i=0,k=nvtxs; i<npes; i++) {
      l = k/(npes-i);
      vtxdist[i+1] = vtxdist[i]+l;
      k -= l;
    }
  }

  MPI_Bcast((void *)vtxdist, npes+1, IDX_DATATYPE, 0, comm);

  graph->gnvtxs = vtxdist[npes];
  graph->nvtxs = vtxdist[mype+1]-vtxdist[mype];
  graph->xadj = idxmalloc(graph->nvtxs+1, "MGridReadTestGraph: xadj");


  if (mype == 0) {
    for (penum=0; penum<npes; penum++) {
      snvtxs = vtxdist[penum+1]-vtxdist[penum];
      sxadj = idxmalloc(snvtxs+1, "MGridReadTestGraph: sxadj");

      idxcopy(snvtxs+1, gxadj+vtxdist[penum], sxadj);
      for (i=snvtxs; i>=0; i--)
        sxadj[i] -= sxadj[0];

      ssize[penum] = gxadj[vtxdist[penum+1]] - gxadj[vtxdist[penum]];

      if (penum == mype) 
        idxcopy(snvtxs+1, sxadj, graph->xadj);
      else
        MPI_Send((void *)sxadj, snvtxs+1, IDX_DATATYPE, penum, 1, comm); 

      free(sxadj);
    }
  }
  else 
    MPI_Recv((void *)graph->xadj, graph->nvtxs+1, IDX_DATATYPE, 0, 1, comm, &status);

  graph->vvol = realmalloc(graph->nvtxs, "MGridReadTestGraph: vvol");

  if (mype == 0) {
    for (penum=0; penum<npes; penum++) {
       snvtxs = vtxdist[penum+1]-vtxdist[penum];

       if (penum == mype)
         realcopy(snvtxs, gvvol+vtxdist[penum], graph->vvol);
       else
         MPI_Send((void *)(gvvol+vtxdist[penum]), snvtxs, REAL_DATATYPE, penum, 1, comm); 
    }
  }
  else 
    MPI_Recv((void *)graph->vvol, graph->nvtxs, REAL_DATATYPE, 0, 1, comm, &status);

  graph->vsurf = realmalloc(graph->nvtxs, "MGridReadTestGraph: vsurf");

  if (mype == 0) {
    for (penum=0; penum<npes; penum++) {
       snvtxs = vtxdist[penum+1]-vtxdist[penum];

       if (penum == mype)
         realcopy(snvtxs, gvsurf+vtxdist[penum], graph->vsurf);
       else
         MPI_Send((void *)(gvsurf+vtxdist[penum]), snvtxs, REAL_DATATYPE, penum, 1, comm); 
    }
  }
  else 
    MPI_Recv((void *)graph->vsurf, graph->nvtxs, REAL_DATATYPE, 0, 1, comm, &status);

  graph->nedges = graph->xadj[graph->nvtxs];
  graph->adjncy = idxmalloc(graph->nedges, "MGridReadTestGraph: graph->adjncy");

  if (mype == 0) {
    for (penum=0; penum<npes; penum++) {

      if (penum == mype) 
        idxcopy(ssize[penum], gadjncy+gxadj[vtxdist[penum]], graph->adjncy);
      else
        MPI_Send((void *)(gadjncy+gxadj[vtxdist[penum]]), ssize[penum], IDX_DATATYPE, penum, 1, comm); 
    }
  }
  else 
    MPI_Recv((void *)graph->adjncy, graph->nedges, IDX_DATATYPE, 0, 1, comm, &status);

  graph->vwgt = NULL;
  graph->adjwgt = realmalloc(graph->nedges, "MGridReadTestGraph: graph->adjwgt");

  if (mype == 0) {
    for (penum=0; penum<npes; penum++) {

      if (penum == mype)
        realcopy(ssize[penum], gadjwgt+gxadj[vtxdist[penum]], graph->adjwgt);
      else
        MPI_Send((void *)(gadjwgt+gxadj[vtxdist[penum]]), ssize[penum], REAL_DATATYPE, penum, 1, comm);
    }

    free(ssize);
  }
  else 
    MPI_Recv((void *)graph->adjwgt, graph->nedges, REAL_DATATYPE, 0, 1, comm, &status);

  if (mype == 0) 
    IMfree(&gxadj, &gvvol, &gvsurf, &gadjncy, &gadjwgt, LTERM);
}


/*************************************************************************
* This function reads the CSR matrix
**************************************************************************/
double *ReadTestCoordinates(MGridGraphType *graph, char *filename, int ndims, MPI_Comm comm)
{
  int i, j, k, npes, mype;
  int penum;
  double *xyz, *txyz;
  FILE *fpin;
  idxtype *vtxdist;
  MPI_Status status;
  char xyzfile[256];

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  vtxdist = graph->vtxdist;

  xyz = fmalloc(graph->nvtxs*ndims, "io");

  if (mype == 0) {
    sprintf(xyzfile, "%s.xyz", filename);
    if ((fpin = fopen(xyzfile, "r")) == NULL) 
      errexit("Failed to open file %s\n", xyzfile);
  }

  if (mype == 0) {
    txyz = fmalloc(2*graph->nvtxs*ndims, "io");

    for (penum=0; penum<npes; penum++) {
      for (k=0, i=vtxdist[penum]; i<vtxdist[penum+1]; i++, k++) {
        for (j=0; j<ndims; j++)
          fscanf(fpin, "%e ", txyz+k*ndims+j);
      }

      if (penum == mype) 
        memcpy((void *)xyz, (void *)txyz, sizeof(float)*ndims*k);
      else {
        MPI_Send((void *)txyz, ndims*k, MPI_FLOAT, penum, 1, comm); 
      }
    }
    free(txyz);
    fclose(fpin);
  }
  else 
    MPI_Recv((void *)xyz, ndims*graph->nvtxs, MPI_FLOAT, 0, 1, comm, &status);

  return xyz;
}



/*************************************************************************
* This function reads the spd matrix
**************************************************************************/
void ReadMGridGraph(char *filename, int *r_nvtxs, idxtype **r_xadj, idxtype **r_adjncy, realtype **r_adjwgt, realtype **r_vvol, realtype **r_vsurf)
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
    nvtxs = 0;
    return;
  }

  sscanf(line, "%d %d", &nvtxs, &nedges);
  nedges *= 2;


  xadj   = idxsmalloc(nvtxs+1, 0, "ReadMGridGraph: xadj");
  vvol   = realmalloc(nvtxs, "ReadMGridGraph: vvol");
  vsurf  = realsmalloc(nvtxs, 0.0, "ReadMGridGraph: vsurf");
  adjncy = idxmalloc(nedges+1, "ReadMGridGraph: adjncy");
  adjwgt = realmalloc(nedges+1, "ReadMGridGraph: adjwgt");


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
    errexit("ReadMGridGraph: Something wrong with the edges from input file %d %d", nedges, k);

  *r_nvtxs = nvtxs;
  *r_xadj = xadj;
  *r_adjncy = adjncy;
  *r_adjwgt = adjwgt;
  *r_vvol = vvol;
  *r_vsurf = vsurf;
}

/*************************************************************************
* This function writes out the partition vector locallly
**************************************************************************/
void WriteParallelPartition(char *fname, idxtype *part, idxtype *vtxdist, int nparts,
                            int mype, int npes)
{
  int i,n;
  char filename[256];
  FILE *fpout;

  n = vtxdist[mype+1] - vtxdist[mype];

  sprintf(filename,"%s.part%d-%d.%d",fname, npes, mype, nparts);

  if ((fpout = fopen(filename, "w")) == NULL)
    errexit("Problems in opening the partition file: %s", filename);

  for (i=0; i<n; i++)
     fprintf(fpout,"%d\n",part[i]);

  fclose(fpout);
}

