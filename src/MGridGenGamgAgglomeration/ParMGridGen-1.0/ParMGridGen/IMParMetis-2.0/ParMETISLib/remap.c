/*
 * remap.c
 *
 * This file contains code that computes the assignment of processors to
 * partition numbers so that it will minimize the redistribution cost
 *
 * George Irene
 */

#include "parmetis.h"

/*************************************************************************
* This function remaps that graph so that it will minimize the 
* redistribution cost
**************************************************************************/
void ReMapGraph(CtrlType *ctrl, GraphType *graph, int rtype, WorkSpaceType *wspace)
{
  int i, nvtxs;
  idxtype *where, *vsize, *map, *lpwgts;

  if (ctrl->npes != ctrl->nparts)
    return;

  IFSET(ctrl->dbglvl, DBG_TIME, MPI_Barrier(ctrl->comm));
  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->RemapTmr));

  nvtxs = graph->nvtxs;
  where = graph->where;
  vsize = graph->vsize;

  map = wspace->pv1;
  lpwgts = idxset(ctrl->npes, 0, wspace->pv2);

  if (vsize == NULL) {
    for (i=0; i<nvtxs; i++)
      lpwgts[where[i]]++;
  }
  else {
    for (i=0; i<nvtxs; i++)
      lpwgts[where[i]] += vsize[i];
  }

  ComputeTotalVReMap1(ctrl, lpwgts, map, wspace);

  for (i=0; i<nvtxs; i++)
    where[i] = map[where[i]];

  IFSET(ctrl->dbglvl, DBG_TIME, MPI_Barrier(ctrl->comm));
  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->RemapTmr));
}


/*************************************************************************
* This function computes the assignment using the the objective the 
* minimization of the total volume of data that needs to move
**************************************************************************/
void ComputeTotalVReMap1(CtrlType *ctrl, idxtype *lpwgts, idxtype *map, WorkSpaceType *wspace)
{
  int i, j, k, npes, mype, nkeys, nmapped, nsaved;
  idxtype mywgt, *rowmap, *orgwgt;
  KeyValueType *order;
  MPI_Status status;

  mype = ctrl->mype;
  npes = ctrl->npes;

  /* Root processor allocates memory for all the lpwgts */
  if (mype == 0) {
    order = (KeyValueType *)IMmalloc(sizeof(KeyValueType)*npes*npes, "ComputeTotalVReMap: order");
    orgwgt = idxmalloc(npes, "ComputeTotalVReMap1: orgwgt");
  }
  else
    order = (KeyValueType *)IMmalloc(sizeof(KeyValueType)*npes, "ComputeTotalVReMap: order");

  for (nkeys=0, i=0; i<npes; i++) {
    if (lpwgts[i] > 0) {
      order[nkeys].key = -lpwgts[i];
      order[nkeys++].val = mype*npes+i;
    }
  }

  /* Gather the weight of the initial assignment */
  mywgt = lpwgts[mype];
  MPI_Gather(&mywgt, 1, IDX_DATATYPE, orgwgt, 1, IDX_DATATYPE, 0, ctrl->comm);

  /* Get into a loop receiving the items in order array */
  if (mype != 0) {
    MPI_Send((void *)order, 2*nkeys, IDX_DATATYPE, 0, 1, ctrl->comm);
  }
  else {
    for (i=1; i<npes; i++) {
      MPI_Recv((void *)(order+nkeys), 2*npes, IDX_DATATYPE, i, 1, ctrl->comm, &status);
      MPI_Get_count(&status, IDX_DATATYPE, &k);
      nkeys += k/2;
    }

    /* Sort them in decreasing order wrt to their absolute values */
    ikeysort(nkeys, order);

    rowmap = idxset(npes, -1, wspace->pv3);
    idxset(npes, -1, map);

    nsaved = 0;
    for (nmapped=k=0; k<nkeys && nmapped<npes; k++) {
      i = order[k].val/npes;
      j = order[k].val%npes;
      if (map[j] == -1 && rowmap[i] == -1) {
        map[j] = i;
        rowmap[i] = j;
        nmapped++;
        nsaved += (-order[k].key - orgwgt[i]);
      }
    }
    
    /* Map unmapped partitions */
    if (nmapped < npes) {
      for (i=j=0; j<npes && nmapped<npes; j++) {
        if (map[j] == -1) {
          for (; i<npes; i++) {
            if (rowmap[i] == -1) {
              map[j] = i;
              rowmap[i] = j;
              nmapped++;
              nsaved += - orgwgt[i];
              break;
            }
          }
        }
      }
    }

    ASSERT(ctrl, nmapped == npes);

    if (nsaved <= 0) {
      for (i=0; i<npes; i++)
        map[i] = i;
    }

    free(orgwgt);
  }

  IMfree(&order, LTERM);

  IFSET(ctrl->dbglvl, DBG_INFO, rprintf(ctrl, "Savings from remapping: %d\n", nsaved)); 

  /* Tell everybody about the map array */
  MPI_Bcast(map, npes, IDX_DATATYPE, 0, ctrl->comm);
 
}
