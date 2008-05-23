/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * mgridgen.h
 *
 * This file includes all necessary header files
 *
 * George Irene
 */


#include <stdio.h>
#ifdef __STDC__
#include <stdlib.h>
#else
/* #include <malloc.h> */
#endif
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <time.h>

#include "IMlib.h"

#ifdef DMALLOC
#include <dmalloc.h>
#else
#include <malloc.h>
#endif

#include "defs.h"
#include "struct.h"
#include "macros.h"
#include "proto.h"
