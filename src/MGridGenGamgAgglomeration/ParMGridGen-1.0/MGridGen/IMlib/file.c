/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * file.c
 *
 * This file contains some simple io functions
 *
 * Irene
 */

#include "IMlib.h"

/*************************************************************************
* This function opens a file
**************************************************************************/
FILE *IMfopen(char *fname, char *mode, char *msg)
{
  FILE *fp;
  char errmsg[256];

  fp = fopen(fname, mode);
  if (fp != NULL)
    return fp;

  sprintf(errmsg,"file: %s, mode: %s, [%s]", fname, mode, msg);
  perror(msg);
  exit(0);
}


/*************************************************************************
* This function closes a file
**************************************************************************/
void IMfclose(FILE *fp)
{
  fclose(fp);
}
