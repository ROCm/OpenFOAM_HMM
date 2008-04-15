/*
 * Copyright 2003, Regents of the University of Minnesota
 *
 * cepic.c
 *
 * This file contains the driving routine for contact/impact simulations
 * for EPIC meshes
 *
 * Started 4/12/03
 * George
 *
 * $Id: cepic-rcb.c,v 1.3 2003/05/03 16:10:48 karypis Exp $
 *
 */

#include <metisbin.h>

#define Flip_int32(type)  (((type >>24) & 0x000000ff) | \
                           ((type >> 8) & 0x0000ff00) | \
                           ((type << 8) & 0x00ff0000) | \
                           ((type <<24) & 0xff000000) )

#define Flip_int64(type)  (((type >>56) & 0x00000000000000ff) | \
                           ((type >>40) & 0x000000000000ff00) | \
                           ((type >>24) & 0x0000000000ff0000) | \
                           ((type >>8)  & 0x00000000ff000000) | \
                           ((type <<8)  & 0x000000ff00000000) | \
                           ((type <<24) & 0x0000ff0000000000) | \
                           ((type <<40) & 0x00ff000000000000) | \
                           ((type <<56) & 0xff00000000000000))


int ComputeMapCost(idxtype nvtxs, idxtype nparts, idxtype *fepart, idxtype *cpart);

/*************************************************************************
* Let the game begin
**************************************************************************/
int main(int argc, char *argv[])
{
  idxtype i, j, istep, options[10], nn, ne, fstep, lstep, nparts, nboxes, u[3], dim, nchanges, ncomm;
  char filename[256];
  idxtype *mien, *mrng, *part, *oldpart, *sflag, *bestdims, *fepart;
  double *mxyz, *bxyz;
  idxtype *xadj, *adjncy, *cntptr, *cntind;
  idxtype numflag = 0, wgtflag = 0, edgecut, etype=2;
  void *cinfo;
  FILE *fpin;
  long long int *ltmp;

  if (argc != 6) {
    mfprintf(stderr, "Usage: %s <nn> <ne> <fstep> <lstep> <nparts>\n", argv[0]);
    exit(0);
  }

  nn     = atoi(argv[1]);
  ne     = atoi(argv[2]);
  fstep  = atoi(argv[3]);
  lstep  = atoi(argv[4]);
  nparts = atoi(argv[5]);

  mprintf("Reading %s, nn: %D, ne: %D, fstep: %D, lstep: %D, nparts: %D\n", filename, nn, ne, fstep, lstep, nparts);

  mien = idxmalloc(4*ne, "main: mien");
  mxyz = gk_dmalloc(3*nn, "main: mxyz");
  mrng = idxmalloc(4*ne, "main: mrng");
  bxyz = gk_dmalloc(6*ne*4, "main: bxyz");

  fepart  = idxmalloc(nn, "main: fepart");
  part    = idxmalloc(nn, "main: part");
  oldpart = idxmalloc(nn, "main: oldpart");
  sflag   = idxmalloc(nn, "main: sflag");

  bestdims  = idxsmalloc(2*nparts, -1, "main: bestdims");

  xadj   = idxmalloc(nn+1, "main: xadj");
  adjncy = idxmalloc(50*nn, "main: adjncy");


  /*========================================================================
   * Read the initial mesh and setup the graph and contact information
   *========================================================================*/
  msprintf(filename, "mien.%04D", fstep);
  fpin = GKfopen(filename, "rb", "main: mien");
  fread(mien, sizeof(int), 4*ne, fpin);
  for (i=0; i<4*ne; i++)
    mien[i] = Flip_int32(mien[i]);
  GKfclose(fpin);

  msprintf(filename, "mxyz.%04D", fstep);
  fpin = GKfopen(filename, "rb", "main: mxyz");
  fread(mxyz, sizeof(double), 3*nn, fpin);
  for (i=0; i<3*nn; i++) {
    ltmp = (long long int *)(mxyz+i);
    *ltmp = Flip_int64(*ltmp);
  }
  GKfclose(fpin);
  mprintf("%e %e %e\n", mxyz[3*0+0], mxyz[3*0+1], mxyz[3*0+2]);

  msprintf(filename, "mrng.%04D", fstep);
  fpin = GKfopen(filename, "rb", "main: mrng");
  fread(mrng, sizeof(int), 4*ne, fpin);
  for (i=0; i<4*ne; i++)
    mrng[i] = Flip_int32(mrng[i]);
  GKfclose(fpin);


  /*========================================================================
   * Determine which nodes are in the surface
   *========================================================================*/
  iset(nn, 0, sflag);
  for (i=0; i<ne; i++) {
    if (mrng[4*i+0] > 0) { /* 1, 2, 3 */
      sflag[mien[4*i+0]-1] = 1;
      sflag[mien[4*i+1]-1] = 1;
      sflag[mien[4*i+2]-1] = 1;
    }
    if (mrng[4*i+1] > 0) { /* 1, 2, 4 */
      sflag[mien[4*i+0]-1] = 1;
      sflag[mien[4*i+1]-1] = 1;
      sflag[mien[4*i+3]-1] = 1;
    }
    if (mrng[4*i+2] > 0) { /* 2, 3, 4 */
      sflag[mien[4*i+1]-1] = 1;
      sflag[mien[4*i+2]-1] = 1;
      sflag[mien[4*i+3]-1] = 1;
    }
    if (mrng[4*i+3] > 0) { /* 1, 3, 4 */
      sflag[mien[4*i+0]-1] = 1;
      sflag[mien[4*i+2]-1] = 1;
      sflag[mien[4*i+3]-1] = 1;
    }
  }

  mprintf("Contact Nodes: %D of %D\n", isum(nn, sflag), nn);


  /*========================================================================
   * Compute the FE partition
   *========================================================================*/
  numflag = mien[idxargmin(4*ne, mien)];
  METIS_MeshToNodal(&ne, &nn, mien, &etype, &numflag, xadj, adjncy);

  options[0] = 0;
  METIS_PartGraphVKway(&nn, xadj, adjncy, NULL, NULL, &wgtflag, &numflag, &nparts,
        options, &edgecut, fepart);

  mprintf("K-way partitioning Volume: %D\n", edgecut);


  /*========================================================================
   * Get into the loop in which you go over the different configurations
   *========================================================================*/
  for (istep=fstep; istep<=lstep; istep++) {
    msprintf(filename, "mxyz.%04D", istep);
    mprintf("Reading %s...............................................................\n", filename);
    fpin = GKfopen(filename, "rb", "main: mxyz");
    fread(mxyz, sizeof(double), 3*nn, fpin);
    for (i=0; i<3*nn; i++) {
      ltmp = (long long int *)(mxyz+i);
      *ltmp = Flip_int64(*ltmp);
    }
    GKfclose(fpin);

    msprintf(filename, "mrng.%04D", istep);
    fpin = GKfopen(filename, "rb", "main: mrng");
    fread(mrng, sizeof(int), 4*ne, fpin);
    for (i=0; i<4*ne; i++)
      mrng[i] = Flip_int32(mrng[i]);
    GKfclose(fpin);

    /* Determine which nodes are in the surface */
    iset(nn, 0, sflag);
    for (i=0; i<ne; i++) {
      if (mrng[4*i+0] > 0) { /* 1, 2, 3 */
        sflag[mien[4*i+0]-1] = 1;
        sflag[mien[4*i+1]-1] = 1;
        sflag[mien[4*i+2]-1] = 1;
      }
      if (mrng[4*i+1] > 0) { /* 1, 2, 4 */
        sflag[mien[4*i+0]-1] = 1;
        sflag[mien[4*i+1]-1] = 1;
        sflag[mien[4*i+3]-1] = 1;
      }
      if (mrng[4*i+2] > 0) { /* 2, 3, 4 */
        sflag[mien[4*i+1]-1] = 1;
        sflag[mien[4*i+2]-1] = 1;
        sflag[mien[4*i+3]-1] = 1;
      }
      if (mrng[4*i+3] > 0) { /* 1, 3, 4 */
        sflag[mien[4*i+0]-1] = 1;
        sflag[mien[4*i+2]-1] = 1;
        sflag[mien[4*i+3]-1] = 1;
      }
    }

    mprintf("Contact Nodes: %D of %D\n", isum(nn, sflag), nn);

    /* Determine the bounding boxes of the surface elements */
    for (nboxes=0, i=0; i<ne; i++) {
      if (mrng[4*i+0] > 0) { /* 1, 2, 3 */
        u[0] = mien[4*i+0]-1;
        u[1] = mien[4*i+1]-1;
        u[2] = mien[4*i+2]-1;
        bxyz[6*nboxes+0] = bxyz[6*nboxes+3] = mxyz[3*u[0]+0];
        bxyz[6*nboxes+1] = bxyz[6*nboxes+4] = mxyz[3*u[0]+1];
        bxyz[6*nboxes+2] = bxyz[6*nboxes+5] = mxyz[3*u[0]+2];
        for (j=1; j<3; j++) {
          for (dim=0; dim<3; dim++) {
            bxyz[6*nboxes+dim] = (bxyz[6*nboxes+dim] > mxyz[3*u[j]+dim] ? mxyz[3*u[j]+dim] : bxyz[6*nboxes+dim]);
            bxyz[6*nboxes+3+dim] = (bxyz[6*nboxes+3+dim] < mxyz[3*u[j]+dim] ? mxyz[3*u[j]+dim] : bxyz[6*nboxes+3+dim]);
          }
        }
        nboxes++;
      }
      if (mrng[4*i+1] > 0) { /* 1, 2, 4 */
        u[0] = mien[4*i+0]-1;
        u[1] = mien[4*i+1]-1;
        u[2] = mien[4*i+3]-1;
        bxyz[6*nboxes+0] = bxyz[6*nboxes+3] = mxyz[3*u[0]+0];
        bxyz[6*nboxes+1] = bxyz[6*nboxes+4] = mxyz[3*u[0]+1];
        bxyz[6*nboxes+2] = bxyz[6*nboxes+5] = mxyz[3*u[0]+2];
        for (j=1; j<3; j++) {
          for (dim=0; dim<3; dim++) {
            bxyz[6*nboxes+dim] = (bxyz[6*nboxes+dim] > mxyz[3*u[j]+dim] ? mxyz[3*u[j]+dim] : bxyz[6*nboxes+dim]);
            bxyz[6*nboxes+3+dim] = (bxyz[6*nboxes+3+dim] < mxyz[3*u[j]+dim] ? mxyz[3*u[j]+dim] : bxyz[6*nboxes+3+dim]);
          }
        }
        nboxes++;
      }
      if (mrng[4*i+2] > 0) { /* 2, 3, 4 */
        u[0] = mien[4*i+1]-1;
        u[1] = mien[4*i+2]-1;
        u[2] = mien[4*i+3]-1;
        bxyz[6*nboxes+0] = bxyz[6*nboxes+3] = mxyz[3*u[0]+0];
        bxyz[6*nboxes+1] = bxyz[6*nboxes+4] = mxyz[3*u[0]+1];
        bxyz[6*nboxes+2] = bxyz[6*nboxes+5] = mxyz[3*u[0]+2];
        for (j=1; j<3; j++) {
          for (dim=0; dim<3; dim++) {
            bxyz[6*nboxes+dim] = (bxyz[6*nboxes+dim] > mxyz[3*u[j]+dim] ? mxyz[3*u[j]+dim] : bxyz[6*nboxes+dim]);
            bxyz[6*nboxes+3+dim] = (bxyz[6*nboxes+3+dim] < mxyz[3*u[j]+dim] ? mxyz[3*u[j]+dim] : bxyz[6*nboxes+3+dim]);
          }
        }
        nboxes++;
      }
      if (mrng[4*i+3] > 0) { /* 1, 3, 4 */
        u[0] = mien[4*i+0]-1;
        u[1] = mien[4*i+2]-1;
        u[2] = mien[4*i+3]-1;
        bxyz[6*nboxes+0] = bxyz[6*nboxes+3] = mxyz[3*u[0]+0];
        bxyz[6*nboxes+1] = bxyz[6*nboxes+4] = mxyz[3*u[0]+1];
        bxyz[6*nboxes+2] = bxyz[6*nboxes+5] = mxyz[3*u[0]+2];
        for (j=1; j<3; j++) {
          for (dim=0; dim<3; dim++) {
            bxyz[6*nboxes+dim] = (bxyz[6*nboxes+dim] > mxyz[3*u[j]+dim] ? mxyz[3*u[j]+dim] : bxyz[6*nboxes+dim]);
            bxyz[6*nboxes+3+dim] = (bxyz[6*nboxes+3+dim] < mxyz[3*u[j]+dim] ? mxyz[3*u[j]+dim] : bxyz[6*nboxes+3+dim]);
          }
        }
        nboxes++;
      }
    }

    cinfo = METIS_PartSurfForContactRCB(&nn, mxyz, sflag, &nparts, part, bestdims);

    METIS_FindContacts(cinfo, &nboxes, bxyz, &nparts, &cntptr, &cntind);

    METIS_FreeContactInfo(cinfo);

    nchanges = 0;
    if (istep > fstep) {
      for (i=0; i<nn; i++)
        nchanges += (part[i] != oldpart[i] ? 1 : 0);
    }
    idxcopy(nn, part, oldpart);

    ncomm = ComputeMapCost(nn, nparts, fepart, part);

    mprintf("Contacting Elements: %D  Indices: %D  Nchanges: %D  MapCost: %D\n", nboxes, cntptr[nboxes]-nboxes, nchanges, ncomm);

    gk_free((void **)&cntptr, &cntind, LTERM);
  }

}  


/***********************************************************************************
* This function determines the cost of moving data between the two meshes assuming
* that a good matching between the two partitions was done!
************************************************************************************/
int ComputeMapCost(idxtype nvtxs, idxtype nparts, idxtype *fepart, idxtype *cpart)
{
  idxtype i, j, k, n, ncomm;
  KeyValueType cand[nparts*nparts];
  idxtype fmatched[nparts], cmatched[nparts];

  /* Compute the overlap */
  for (i=0; i<nparts; i++) {
    for (j=0; j<nparts; j++) {
      cand[i*nparts+j].key = 0;
      cand[i*nparts+j].val = i*nparts+j;
    }
  }

  for (k=0, i=0; i<nvtxs; i++) {
    if (cpart[i] >= 0) {
      cand[(fepart[i]-1)*nparts+(cpart[i]-1)].key++;
      k++;
    }
  }

mprintf("Contact points: %D\n", k);
      
  ikeysort(nparts*nparts, cand);

  iset(nparts, -1, fmatched);
  iset(nparts, -1, cmatched);


  for (ncomm=0, k=nparts*nparts-1; k>=0; k--) {
    i = cand[k].val/nparts;
    j = cand[k].val%nparts;

    if (fmatched[i] == -1 && cmatched[j] == -1) {
      fmatched[i] = j;
      cmatched[j] = i;
    }
    else 
      ncomm += cand[k].key;
  }

mprintf("Ncomm: %D\n", ncomm);

  return ncomm;

}


