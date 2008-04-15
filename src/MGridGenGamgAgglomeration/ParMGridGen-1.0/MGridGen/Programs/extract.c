/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * extract.c
 *
 * This file is the driver routine for extracting the graph
 * for each of the fused elements I used these graphs to call onmetis
 *
 * Irene
 */

#include "mgridgen.h"

void WriteGraph(int, idxtype *, idxtype *, idxtype, realtype, int, int);
void  ExtractGraph(int, idxtype *, idxtype *, idxtype *, realtype *, realtype *,
                   realtype *, int, idxtype *);

/*************************************************************************
* Let the game begin
**************************************************************************/
int main(int argc, char *argv[])
{
  int nvtxs, nedges;
  int options[10], nmoves, nparts, minsize, maxsize;
  idxtype *part;
  idxtype *xadj, *adjncy, *vwgt;
  realtype *vvol, *vsurf, *adjwgt;
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
    printf("\tminsize:\tA lower bound on the cell size (suggested)\n");
    printf("\tmaxsize:\tAn upper bound on the cell size (strict)\n");
    printf("--------------------------------------------------------------------\n");
    printf("Recomended usage: %s <GraphFile> <Dim> 4 6 ? ? 32\n", argv[0]);
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

  printf("Parameters: %s, Dim=%d [%2d %2d] CType=%d RType=%d Nvtxs=%d Nedges=%d\n", filename,
          options[OPTION_DIM], minsize, maxsize, options[OPTION_CTYPE], options[OPTION_RTYPE],
          graph.nvtxs, graph.xadj[graph.nvtxs]);

  /* Keep graph info */
  nvtxs = graph.nvtxs;
  nedges = graph.xadj[nvtxs];
  part = imalloc(nvtxs, "main: part");
  xadj = imalloc(nvtxs+1, "main: xadj");
  adjncy = imalloc(nedges, "main: adjncy");
  vwgt = imalloc(nvtxs, "main: vwgt");
  vvol = fmalloc(nvtxs, "main:vvol");
  vsurf = fmalloc(nvtxs, "main:vsurf");
  adjwgt = fmalloc(nedges, "main:adjwgt");
  icopy(nvtxs+1, graph.xadj, xadj);
  icopy(nedges, graph.adjncy, adjncy);
  iset(nvtxs, 1, vwgt);
  fcopy(nvtxs, graph.vvol, vvol);
  fcopy(nvtxs, graph.vsurf, vsurf);
  fcopy(nedges, graph.adjwgt, adjwgt);

  cleartimer(tmr);
  starttimer(tmr);
  MGridGen(graph.nvtxs, graph.xadj, graph.vvol, graph.vsurf, graph.adjncy, graph.adjwgt, 
           minsize, maxsize, options, &nmoves, &nparts, part);
  stoptimer(tmr);
  printf("Total Time: %lf\n", gettimer(tmr));

  WritePartition("partitions", part, graph.nvtxs, nparts);

  ExtractGraph(nvtxs, xadj, adjncy, vwgt, vvol, vsurf, adjwgt, nparts, part);

  IMfree(&graph.xadj, &graph.vvol, &graph.vsurf, &graph.adjncy, &graph.adjwgt, &part, LTERM);
  return(0);
}  

/***************************************************************************/
void  ExtractGraph(int nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
                   realtype *vvol, realtype *vsurf, realtype *adjwgt,
                   int nparts, idxtype *part)
{

   int i, j, k, l, last, ipart;
   int nfvtxs, nedges;
   idxtype *map, *fxadj, *fadjncy;
   idxtype *ptr, *ind, *pwgts;
   realtype aratio;
   realtype *pvol, *psurf;

  /* Create a part-to-vertex mapping; compute pwgts, psurf, pvol */
  ptr   = ismalloc(nparts+1, 0, "ExtractGraph: ptr");
  ind   = imalloc(nvtxs, "ExtractGraph: ind");
  pwgts = ismalloc(nparts, 0, "ExtractGraph: pwgts");
  pvol  = fsmalloc(nparts, 0.0, "ExtractGraph: pvol");
  psurf = fsmalloc(nparts, 0.0, "ExtractGraph: psurf");

  for (i=0; i<nvtxs; i++) {
     ipart = part[i];
     ptr[ipart]++;
     pwgts[ipart] += vwgt[i];
     pvol[ipart] += vvol[i];
     psurf[ipart] += vsurf[i];
     for (j=xadj[i]; j<xadj[i+1]; j++) {
        if (part[adjncy[j]] != ipart)
          psurf[ipart] += adjwgt[j];
     }
  }

  nfvtxs = iamax(nparts, ptr);
  nfvtxs = ptr[nfvtxs];
  nedges = xadj[nvtxs];

  MAKECSR(i, nparts, ptr);
  for (i=0; i<nvtxs; i++) {
     ipart = part[i];
     ind[ptr[ipart]] = i;
     ptr[ipart]++;
  }

  for (ipart=nparts; ipart>0; ipart--)
     ptr[ipart] = ptr[ipart-1];
  ptr[0] = 0;


  /* Create the connectivity graph for each of the fused elements */
  fxadj = ismalloc(nfvtxs+1, 0, "FusedElementGraph: fxadj");
  fadjncy = imalloc(nedges, "FusedElementGraph: fadjncy");
  map = ismalloc(nvtxs, -1, "FusedElementGraph: map"); 

  for (ipart=0; ipart<nparts; ipart++) {
     for (nfvtxs=0, l=ptr[ipart]; l<ptr[ipart+1]; l++, nfvtxs++) 
        map[ind[l]] = nfvtxs;

     for (fxadj[0]=last=0, l=ptr[ipart]; l<ptr[ipart+1]; l++) {
        i = ind[l];
        for (j=xadj[i]; j<xadj[i+1]; j++) {
           k=adjncy[j];
           if (part[k] == ipart && k != i ) 
             fadjncy[last++] = map[k];
        }
        fxadj[l+1-ptr[ipart]] = last;
     }

     aratio = ARATIO(3, psurf[ipart], pvol[ipart]);
     WriteGraph(nfvtxs, fxadj, fadjncy, pwgts[ipart], aratio, nparts, ipart);

     iset(nvtxs, -1, map);
  }

  IMfree(ptr, ind, fxadj, fadjncy, map, LTERM);
}

void WriteGraph(int nvtxs, idxtype *xadj, idxtype *adjncy, idxtype pwgts,
                realtype aratio, int nparts, int mypart)
{
   int i, j;
   char filename[20];
   FILE *fp;

   sprintf(filename,"test.graph%d.%d",nparts,mypart);
   fp = fopen(filename, "w");

   fprintf(fp, "%d %d", nvtxs, xadj[nvtxs]/2);

   for (i=0; i<nvtxs; i++) {
      fprintf(fp, "\n");
      for (j=xadj[i]; j<xadj[i+1]; j++)
        fprintf(fp, " %d", adjncy[j]+1);
  }

  fprintf(fp, "\n pwgts=%d aratio=%f\n", pwgts, aratio);

  fclose(fp);
}
