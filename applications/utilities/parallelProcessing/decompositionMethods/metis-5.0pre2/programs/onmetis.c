/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * onmetis.c
 *
 * This file contains the driving routine for multilevel method
 *
 * Started 8/28/94
 * George
 *
 * $Id: onmetis.c,v 1.2 2002/08/10 06:02:53 karypis Exp $
 *
 */

#include <metisbin.h>



/*************************************************************************
* Let the game begin
**************************************************************************/
int main(int argc, char *argv[])
{
  idxtype i, options[10];
  idxtype *perm, *iperm;
  GraphType graph;
  char filename[256];
  idxtype numflag = 0, wgtflag;
  double TOTALTmr, METISTmr, IOTmr, SMBTmr;


  if (argc != 2) {
    mprintf("Usage: %s <GraphFile>\n",argv[0]);
    exit(0);
  }
    
  strcpy(filename, argv[1]);

  gk_clearcputimer(TOTALTmr);
  gk_clearcputimer(METISTmr);
  gk_clearcputimer(IOTmr);
  gk_clearcputimer(SMBTmr);

  gk_startcputimer(TOTALTmr);
  gk_startcputimer(IOTmr);
  ReadGraph(&graph, filename, &wgtflag);
  if (graph.nvtxs <= 0) {
    mprintf("Empty graph. Nothing to do.\n");
    exit(0);
  }
  if (graph.ncon != 1) {
    mprintf("Ordering can only be applied to graphs with one constraint.\n");
    exit(0);
  }
  gk_stopcputimer(IOTmr);

  /* Ordering does not use weights! */
  gk_free((void **)&graph.vwgt, &graph.adjwgt, LTERM);

  mprintf("**********************************************************************\n");
  mprintf("%s", METISTITLE);
  mprintf("Graph Information ---------------------------------------------------\n");
  mprintf("  Name: %s, #Vertices: %D, #Edges: %D\n\n", filename, graph.nvtxs, graph.nedges/2);
  mprintf("Node-Based Ordering... ----------------------------------------------\n");

  perm = idxmalloc(graph.nvtxs, "main: perm");
  iperm = idxmalloc(graph.nvtxs, "main: iperm");
  options[0] = 0;

  gk_startcputimer(METISTmr);
  METIS_NodeND(&graph.nvtxs, graph.xadj, graph.adjncy, &numflag, options, perm, iperm);
  gk_stopcputimer(METISTmr);

  gk_startcputimer(IOTmr);
  WritePermutation(filename, iperm, graph.nvtxs); 
  gk_stopcputimer(IOTmr);

  gk_startcputimer(SMBTmr);
  ComputeFillIn(&graph, iperm);
  gk_stopcputimer(SMBTmr);

  gk_stopcputimer(TOTALTmr);

  mprintf("\nTiming Information --------------------------------------------------\n");
  mprintf("  I/O:                     \t %7.3f\n", gk_getcputimer(IOTmr));
  mprintf("  Ordering:                \t %7.3f   (ONMETIS time)\n", gk_getcputimer(METISTmr));
  mprintf("  Symbolic Factorization:  \t %7.3f\n", gk_getcputimer(SMBTmr));
  mprintf("  Total:                   \t %7.3f\n", gk_getcputimer(TOTALTmr));
  mprintf("**********************************************************************\n");


  gk_free((void **)&graph.xadj, &graph.adjncy, &perm, &iperm, LTERM);
}  


