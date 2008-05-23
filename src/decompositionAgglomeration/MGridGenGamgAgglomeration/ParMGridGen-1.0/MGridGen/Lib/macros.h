/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * macros.h
 *
 * This file contains macros used in multilevel
 *
 * George Irene
 */

#define ARATIO1(dim, surf, vol) ((dim == 2) ? (pow((surf), 2)/(vol)) : (pow((surf), 1.5)/(vol)))
#define ARATIO(dim, surf, vol)  ((dim == 2) ? ((surf)*(surf)/(vol)) : (sqrt((surf)*(surf)*(surf))/(vol)))
#define ARATIO2(dim, surf, vol) ((dim == 2) ? ((surf)*(surf)*(surf)*(surf)/(vol)*(vol)) : ((surf)*(surf)*(surf)/((vol)*(vol))))

/*************************************************************************
* This macro is used to handle dbglvl
**************************************************************************/
#define IFSET(a, flag, cmd) if ((a)&(flag)) (cmd);
