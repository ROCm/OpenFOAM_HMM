/*
 * struct.h
 *
 * This file contains data structures for the various programs of METIS.
 *
 * Started 8/9/02
 * George
 *
 * $Id: struct.h,v 1.8 2003/04/04 23:22:50 karypis Exp $
 */

#ifndef _STRUCTBIN_H_
#define _STRUCTBIN_H_


/*************************************************************************
* The following data structure implements a string-2-idxtype mapping
* table used for parsing command-line options
**************************************************************************/
typedef struct {
  char *name;
  idxtype id;
} StringMapType;



/*************************************************************************
* This data structure stores the various command line arguments
**************************************************************************/
typedef struct {
  idxtype mtype;
  idxtype itype;
  idxtype rtype;

  idxtype balance;
  idxtype ntrials;
  idxtype niter;

  idxtype seed;
  idxtype dbglvl;

  idxtype nparts;

  char *filename;
  char *xyzfilename;
  char *tpwgts;

  float iotimer;
  float clustertimer;
  float reporttimer;
} ParamType;


#endif 
