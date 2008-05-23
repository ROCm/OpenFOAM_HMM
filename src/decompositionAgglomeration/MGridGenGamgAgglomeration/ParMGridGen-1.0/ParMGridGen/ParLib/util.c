/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * util.c
 *
 * This function contains various utility routines
 *
 * George Irene
 */

#include "parmgridgen.h"


/*************************************************************************
* This function prints an error message and exits
**************************************************************************/
void MGriderrexit(char *f_str,...)
{
  va_list argp;
  char out1[256], out2[256];

  va_start(argp, f_str);
  vsprintf(out1, f_str, argp);
  va_end(argp);

  sprintf(out2, "Error! %s", out1);

  fprintf(stdout, out2);
  fflush(stdout);

  abort();
}


/*************************************************************************
* This function prints an error message and exits
**************************************************************************/
void MGridmyprintf(MGridCtrlType *ctrl, char *f_str,...)
{
  va_list argp;
  char out1[256], out2[256];

  va_start(argp, f_str);
  vsprintf(out1, f_str, argp);
  va_end(argp);

  sprintf(out2, "[%2d] %s", ctrl->mype, out1);

  fprintf(stdout, out2);
  fflush(stdout);

}



/*************************************************************************
* This function prints an error message and exits
**************************************************************************/
void MGridrprintf(MGridCtrlType *ctrl, char *f_str,...)
{
  va_list argp;

  if (ctrl->mype == 0) {
    va_start(argp, f_str);
    vfprintf(stdout, f_str, argp);
    va_end(argp);
  }

  fflush(stdout);

  MPI_Barrier(ctrl->comm);

}

#ifdef XXXX
/*************************************************************************
* This function returns the seconds
**************************************************************************/
float seconds(void)
{
  return ((float) clock()/CLOCKS_PER_SEC);
}
#endif
