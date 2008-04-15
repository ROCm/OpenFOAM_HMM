/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * graphchk.c
 *
 * This file checks the validity of a graph
 *
 * Started 8/28/94
 * George
 *
 * $Id: graphchk.c,v 1.2 2002/08/10 06:02:53 karypis Exp $
 *
 */

#include <metisbin.h>



/*************************************************************************
* Let the game begin
**************************************************************************/
int main(int argc, char *argv[])
{
  GraphType graph;
  char filename[256];
  idxtype wgtflag;

  if (argc != 2) {
    mprintf("Usage: %s <GraphFile>\n", argv[0]);
    exit(0);
  }
    
  strcpy(filename, argv[1]);

  ReadGraph(&graph, filename, &wgtflag);
  if (graph.nvtxs == 0) {
    mprintf("Empty graph!\n");
    exit(0);
  }

  mprintf("**********************************************************************\n");
  mprintf("%s", METISTITLE);
  mprintf("Graph Information ---------------------------------------------------\n");
  mprintf("  Name: %s, #Vertices: %D, #Edges: %D\n\n", filename, graph.nvtxs, graph.nedges/2);
  mprintf("Checking Graph... ---------------------------------------------------\n");

  if (CheckGraph(&graph))
    mprintf("   The format of the graph is correct!\n");
  else
    mprintf("   The format of the graph is incorrect!\n");

  mprintf("\n**********************************************************************\n");


  gk_free((void **)&graph.xadj, &graph.adjncy, &graph.vwgt, &graph.adjwgt, LTERM);
}  


