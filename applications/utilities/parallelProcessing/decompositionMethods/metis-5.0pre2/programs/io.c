/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * io.c
 *
 * This file contains routines related to I/O
 *
 * Started 8/28/94
 * George
 *
 * $Id: io.c,v 1.7 2003/04/13 04:45:14 karypis Exp $
 *
 */

#include <metisbin.h>



/*************************************************************************
* This function reads the spd matrix
**************************************************************************/
void ReadGraph(GraphType *graph, char *filename, idxtype *wgtflag)
{
  idxtype i, j, k, l, fmt, readew, readvw, ncon, edge, ewgt;
  idxtype *xadj, *adjncy, *vwgt, *adjwgt;
  char *line, *oldstr, *newstr;
  FILE *fpin;

  InitGraph(graph);

  line = gk_cmalloc(MAXLINE+1, "ReadGraph: line");

  fpin = gk_fopen(filename, "r", __func__);

  do {
    fgets(line, MAXLINE, fpin);
  } while (line[0] == '%' && !feof(fpin));

  if (feof(fpin)) {
    graph->nvtxs = 0;
    gk_free((void **)&line, LTERM);
    return;
  }

  fmt = ncon = 0;
  msscanf(line, "%D %D %D %D", &(graph->nvtxs), &(graph->nedges), &fmt, &ncon);


  readew = (fmt%10 > 0);
  readvw = ((fmt/10)%10 > 0);
  if (fmt >= 100) {
    mprintf("Cannot read this type of file format!");
    exit(0);
  }


  *wgtflag = 0;
  if (readew)
    *wgtflag += 1;
  if (readvw)
    *wgtflag += 2;

  if (ncon > 0 && !readvw) {
    mprintf("------------------------------------------------------------------------------\n");
    mprintf("***  I detected an error in your input file  ***\n\n");
    mprintf("You specified ncon=%D, but the fmt parameter does not specify vertex weights\n", ncon);
    mprintf("Make sure that the fmt parameter is set to either 10 or 11.\n");
    mprintf("------------------------------------------------------------------------------\n");
    exit(0);
  }

  graph->nedges *=2;
  ncon = graph->ncon = (ncon == 0 ? 1 : ncon);

  /*mprintf("%D %D %D %D %D [%D %D]\n", fmt, fmt%10, (fmt/10)%10, ncon, graph->ncon, readew, readvw);*/

  if (graph->nvtxs > MAXIDX) 
    errexit("\nThe matrix is too big: %d [%d %d]\n", graph->nvtxs, MAXIDX, sizeof(idxtype));

  xadj = graph->xadj = idxsmalloc(graph->nvtxs+1, 0, "ReadGraph: xadj");
  adjncy = graph->adjncy = idxmalloc(graph->nedges, "ReadGraph: adjncy");

  vwgt = graph->vwgt = (readvw ? idxmalloc(ncon*graph->nvtxs, "ReadGraph: vwgt") : NULL);
  adjwgt = graph->adjwgt = (readew ? idxmalloc(graph->nedges, "ReadGraph: adjwgt") : NULL);

  /* Start reading the graph file */
  for (xadj[0]=0, k=0, i=0; i<graph->nvtxs; i++) {
    do {
      fgets(line, MAXLINE, fpin);
    } while (line[0] == '%' && !feof(fpin));
    oldstr = line;
    newstr = NULL;

    if (strlen(line) == MAXLINE) 
      errexit("\nBuffer for fgets not big enough!\n");

    if (readvw) {
      for (l=0; l<ncon; l++) {
        vwgt[i*ncon+l] = strtoidx(oldstr, &newstr, 10);
        oldstr = newstr;
      }
    }

    for (;;) {
      edge = strtoidx(oldstr, &newstr, 10) -1;
      oldstr = newstr;

      if (readew) {
        ewgt = strtoidx(oldstr, &newstr, 10);
        oldstr = newstr;
      }

      if (edge < 0)
        break;

      adjncy[k] = edge;
      if (readew) 
        adjwgt[k] = ewgt;
      k++;
    } 
    xadj[i+1] = k;
  }

  gk_fclose(fpin);

  if (k != graph->nedges) {
    mprintf("------------------------------------------------------------------------------\n");
    mprintf("***  I detected an error in your input file  ***\n\n");
    mprintf("In the first line of the file, you specified that the graph contained\n%D edges. However, I only found %D edges in the file.\n", graph->nedges/2, k/2);
    if (2*k == graph->nedges) {
      mprintf("\n *> I detected that you specified twice the number of edges that you have in\n");
      mprintf("    the file. Remember that the number of edges specified in the first line\n");
      mprintf("    counts each edge between vertices v and u only once.\n\n");
    }
    mprintf("Please specify the correct number of edges in the first line of the file.\n");
    mprintf("------------------------------------------------------------------------------\n");
    exit(0);
  }

  gk_free((void **)&line, LTERM);
}


/*************************************************************************
* This function reads the spd matrix
**************************************************************************/
void ReadCoordinates(GraphType *graph, char *filename)
{
  idxtype i, j, k, l, nvtxs, fmt, readew, readvw, ncon, edge, ewgt;
  FILE *fpin;
  char *line;


  fpin = gk_fopen(filename, "r", __func__);

  nvtxs = graph->nvtxs;

  graph->coords = gk_dsmalloc(3*nvtxs, 0.0, "ReadCoordinates: coords");

  line = gk_cmalloc(MAXLINE+1, "ReadCoordinates: line");

  for (i=0; i<nvtxs; i++) {
    fgets(line, MAXLINE, fpin);
    msscanf(line, "%lf %lf %lf", graph->coords+3*i+0, graph->coords+3*i+1, graph->coords+3*i+2);
  }
    
  gk_fclose(fpin);

  gk_free((void **)&line, LTERM);
}


/*************************************************************************
* This function writes out the partition vector
**************************************************************************/
void WritePartition(char *fname, idxtype *part, idxtype n, idxtype nparts)
{
  FILE *fpout;
  idxtype i;
  char filename[256];

  msprintf(filename,"%s.part.%D",fname, nparts);

  fpout = gk_fopen(filename, "w", __func__);

  for (i=0; i<n; i++)
    fprintf(fpout,"%" PRIIDX "\n", part[i]);

  gk_fclose(fpout);

}


/*************************************************************************
* This function writes out the partition vectors for a mesh
**************************************************************************/
void WriteMeshPartition(char *fname, idxtype nparts, idxtype ne, idxtype *epart, 
       idxtype nn, idxtype *npart)
{
  FILE *fpout;
  idxtype i;
  char filename[256];

  msprintf(filename,"%s.epart.%D",fname, nparts);

  fpout = gk_fopen(filename, "w", __func__);

  for (i=0; i<ne; i++)
    fprintf(fpout,"%" PRIIDX "\n", epart[i]);

  gk_fclose(fpout);


  msprintf(filename,"%s.npart.%D",fname, nparts);

  fpout = gk_fopen(filename, "w", __func__);

  for (i=0; i<nn; i++)
    fprintf(fpout, "%" PRIIDX "\n", npart[i]);

  gk_fclose(fpout);

}



/*************************************************************************
* This function writes out the partition vector
**************************************************************************/
void WritePermutation(char *fname, idxtype *iperm, idxtype n)
{
  FILE *fpout;
  idxtype i;
  char filename[256];

  msprintf(filename,"%s.iperm",fname);

  fpout = gk_fopen(filename, "w", __func__);

  for (i=0; i<n; i++)
    fprintf(fpout, "%" PRIIDX "\n", iperm[i]);

  gk_fclose(fpout);

}



/*************************************************************************
* This function checks if a graph is valid
**************************************************************************/
int CheckGraph(GraphType *graph)
{
  idxtype i, j, k, l, nvtxs, err=0;
  idxtype *xadj, *adjncy, *adjwgt;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;


  for (i=0; i<nvtxs; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];

      if (i == k) {
        mprintf("Vertex %D contains a self-loop (i.e., diagonal entry in the matrix)!\n", i);
        err++;
      }
      else {
        for (l=xadj[k]; l<xadj[k+1]; l++) {
          if (adjncy[l] == i) {
            if (adjwgt != NULL && adjwgt[l] != adjwgt[j]) {
              mprintf("Edges (%D %D) and (%D %D) do not have the same weight! %D %D\n", i,k,k,i, adjwgt[l], adjwgt[adjncy[j]]);
              err++;
            }
            break;
          }
        }
        if (l == xadj[k+1]) {
          mprintf("Missing edge: (%D %D)!\n", k, i);
          err++;
        }
      }
    }
  }

  if (err > 0) 
    mprintf("A total of %D errors exist in the input file. Correct them, and run again!\n", err);

  return (err == 0 ? 1 : 0);
}


/****************************************************************************
* This function detect the input mesh type
***************************************************************************/
int MeshType(char *filename)
{
  int i, j, k, l, len, cnt=0;
  FILE *fpin;
  char temp[40], inpt[80];
  int firstline[3];

  fpin = gk_fopen(filename, "r", __func__);

  mfscanf(fpin,"%[^\n]s", inpt);

  gk_fclose(fpin); 

  len = strlen(inpt);

  i=0;k=0;
  while (inpt[i]==' ') i++;
  while (i<=len) {
    if (inpt[i]==' ' || i==len) {
      l=0;
      for (j=k; j<i;j++ )
        temp[l++]=inpt[j];

      temp[l]='\0';
      firstline[cnt++] = atoi(temp);
      while (inpt[i]==' ') i++;
      k=i;
      if (i==len) break;
    }
    else 
      i++;
  }


  if (cnt==1) 
    return 0;     /*Mixed element without weight */
  else if (cnt==2 && firstline[1]>0) 
    return 1;  /*Fixed element without weight*/
  else if (cnt==2 && firstline[1]==-1) 
    return 2; /*Mixed element with weight*/ 
  else if (cnt==3 && firstline[2]==-1) 
    return 3; /*Fixed element with weight*/ 

}



/*************************************************************************
* This function reads the element node array of a mesh
**************************************************************************/
idxtype *ReadMesh(char *filename, idxtype *ne, idxtype *nn, idxtype *etype)
{
  idxtype i, j, k, esize;
  idxtype *elmnts;
  FILE *fpin;

  fpin = gk_fopen(filename, "r", __func__);

  mfscanf(fpin, "%D %D", ne, etype);

  switch (*etype) {
    case 1:
      esize = 3;
      break;
    case 2:
      esize = 4;
      break;
    case 3:
      esize = 8;
      break;
    case 4:
      esize = 4;
      break;
    case 5:
      esize = 2;
      break;
    default:
      errexit("Unknown mesh-element type: %d\n", *etype);
  }

  elmnts = idxmalloc(esize*(*ne), "ReadMesh: elmnts");

  for (j=esize*(*ne), i=0; i<j; i++) {
    mfscanf(fpin, "%D", elmnts+i);
    elmnts[i]--;
  }

  gk_fclose(fpin);

  *nn = elmnts[idxargmax(j, elmnts)]+1;

  return elmnts;
}


/*************************************************************************
* This function reads the element node array of a mesh with weight
**************************************************************************/
idxtype *ReadMeshWgt(char *filename, idxtype *ne, idxtype *nn, idxtype *etype, 
            idxtype *vwgt)
{
  idxtype i, j, k, esize;
  idxtype *elmnts;
  FILE *fpin;

  fpin = gk_fopen(filename, "r", __func__);

  mfscanf(fpin, "%D %D", ne, etype);
  mfscanf(fpin, "%D", &i);

  switch (*etype) {
    case 1:
      esize = 3;
      break;
    case 2:
      esize = 4;
      break;
    case 3:
      esize = 8;
      break;
    case 4:
      esize = 4;
      break;
    case 5:
      esize = 2;
      break;
    default:
      errexit("Unknown mesh-element type: %d\n", *etype);
  }

  elmnts = idxmalloc(esize*(*ne), "ReadMeshWgt: elmnts");

  for (j=0, i=0; i<*ne; i++) {
    mfscanf(fpin, "%D", vwgt+i);
    for (k=0; k<esize; k++) {
      mfscanf(fpin, "%D", elmnts+j);
      elmnts[j++]--;
    }
 }

  gk_fclose(fpin);

  *nn = elmnts[idxargmax(j, elmnts)]+1;

  return elmnts;
}


/*************************************************************************
* This function reads the weights of each elements
**************************************************************************/
idxtype *ReadWgt(char *filename, idxtype *ne, idxtype *nn, idxtype *etype)
{
  idxtype i, j, k, l, esize;
  idxtype *vwgt;
  FILE *fpin;

  fpin = gk_fopen(filename, "r", __func__);

  mfscanf(fpin, "%D %D", ne, etype);
  mfscanf(fpin, "%D", &i);

  switch (*etype) {
    case 1:
      esize = 3;
      break;
    case 2:
      esize = 4;
      break;
    case 3:
      esize = 8;
      break;
    case 4:
      esize = 4;
      break;
    case 5:
      esize = 2;
      break;
    default:
      errexit("Unknown mesh-element type: %d\n", *etype);
  }

  vwgt = idxmalloc(*ne, "ReadWgt: vwgt");

  for (j=0, i=0; i<*ne; i++) {
    mfscanf(fpin, "%D", vwgt+i);
    for (k=0; k<esize; k++) {
      mfscanf(fpin, "%D", &l);
      j++;
    }
  }

  gk_fclose(fpin);

  return vwgt;
}


/*************************************************************************
* This function reads # of  element of a mixed mesh
**************************************************************************/
idxtype MixedElements(char *filename)
{
  idxtype ne;
  FILE *fpin;

  fpin = gk_fopen(filename, "r", __func__);

  mfscanf(fpin, "%D", &ne);
  gk_fclose(fpin);

  return ne;
}

/*************************************************************************
* This function reads the element node array of a i Mixed mesh
**************************************************************************/
idxtype *ReadMixedMesh(char *filename, idxtype *ne, idxtype *nn, idxtype  *etype)
{
  idxtype i, j, k, esize;
  idxtype *elmnts;
  FILE *fpin;
  idxtype sizes[]={-1,3,4,8,4,2};

  fpin = gk_fopen(filename, "r", __func__);

  mfscanf(fpin, "%D", ne);

  elmnts = idxmalloc(8*(*ne), "ReadMixedMesh: elmnts");
  
  for (j=0, i=0; i<*ne; i++) {
    mfscanf(fpin, "%D", etype+i);
    for (k=0;k<sizes[etype[i]];k++) {
      mfscanf(fpin, "%D", elmnts+j);
      elmnts[j++]--;
    }
  }

  gk_fclose(fpin);

  *nn = elmnts[idxargmax(j, elmnts)]+1;

  return elmnts;
}


/*************************************************************************
* This function reads the element node array of a  Mixed mesh with weight
**************************************************************************/
idxtype *ReadMixedMeshWgt(char *filename, idxtype *ne, idxtype *nn, 
           idxtype  *etype, idxtype *vwgt)
{
  idxtype i, j, k, esize;
  idxtype *elmnts;
  FILE *fpin;
  idxtype sizes[]={-1,3,4,8,4,2};

  fpin = gk_fopen(filename, "r", __func__);

  mfscanf(fpin, "%D", ne);

  mfscanf(fpin, "%D", nn);

  elmnts = idxmalloc(8*(*ne), "ReadMixedMeshWgt: elmnts");

  for (j=0, i=0; i<*ne; i++) {
    mfscanf(fpin, "%D",etype+i);
    mfscanf(fpin, "%D",vwgt+i);

    for (k=0;k<sizes[etype[i]];k++) {
      mfscanf(fpin, "%D", elmnts+j);
      elmnts[j++]--;
    }
  }

  gk_fclose(fpin);

  *nn = elmnts[idxargmax(j, elmnts)]+1;

  return elmnts;
}

/************************************************************************
* This function reads the element node array of a i Mixed mesh
**************************************************************************/
idxtype *ReadMgcnums(char *filename)
{
  idxtype i;
  idxtype *mgc;
  FILE *fpin;

  fpin = gk_fopen(filename, "r", __func__);

  mgc = idxmalloc(36, "Readmgcnums: mgcnums");

  for (i=0; i<36; i++) {
    if (i<6 || i%6==0) 
      mgc[i]=-1;
    else 
      mfscanf(fpin, "%D", mgc+i);
  }

  gk_fclose(fpin);

  return mgc;
}






/*************************************************************************
* This function writes a graphs into a file 
**************************************************************************/
void WriteGraph(char *filename, idxtype nvtxs, idxtype *xadj, idxtype *adjncy)
{
  idxtype i, j;
  FILE *fpout;

  fpout = gk_fopen(filename, "w", __func__);

  mfprintf(fpout, "%D %D", nvtxs, xadj[nvtxs]/2);

  for (i=0; i<nvtxs; i++) {
    mfprintf(fpout, "\n");
    for (j=xadj[i]; j<xadj[i+1]; j++)
      fprintf(fpout, " %" PRIIDX, adjncy[j]+1);
  }

  gk_fclose(fpout);
}


/*************************************************************************
* This function writes weighted  graph into a file
**************************************************************************/
void WriteWgtGraph(char *filename, idxtype nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt)
{
  idxtype i, j;
  FILE *fpout;

  fpout = gk_fopen(filename, "w", __func__);

  mfprintf(fpout, "%D %D", nvtxs, xadj[nvtxs]/2);
  mfprintf(fpout, " %D", 10);

  for (i=0; i<nvtxs; i++) {
    fprintf(fpout, "\n");
    fprintf(fpout, "%" PRIIDX, vwgt[i]);
    
    for (j=xadj[i]; j<xadj[i+1]; j++)
      fprintf(fpout, " %" PRIIDX, adjncy[j]+1);
  }

  gk_fclose(fpout);
}

/*************************************************************************
* This function writes a graphs into a file 
**************************************************************************/
void WriteMocGraph(GraphType *graph)
{
  idxtype i, j, nvtxs, ncon;
  idxtype *xadj, *adjncy;
  float *nvwgt;
  char filename[256];
  FILE *fpout;

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  nvwgt = graph->nvwgt;

  msprintf(filename, "moc.graph.%D.%D", nvtxs, ncon);

  fpout = gk_fopen(filename, "w", __func__);

  mfprintf(fpout, "%D %D 10 1 %D", nvtxs, xadj[nvtxs]/2, ncon);

  for (i=0; i<nvtxs; i++) {
    mfprintf(fpout, "\n");
    for (j=0; j<ncon; j++)
      fprintf(fpout, "%" PRIIDX " ", (int)((float)10e6*nvwgt[i*ncon+j]));

    for (j=xadj[i]; j<xadj[i+1]; j++)
      fprintf(fpout, " %" PRIIDX, adjncy[j]+1);
  }

  gk_fclose(fpout);
}
