/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * parmgridgen.h
 *
 * This file includes all necessary header files
 *
 * George Irene
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <time.h>
#include <mpi.h>

#ifdef DMALLOC
#include <dmalloc.h>
#else
#include <malloc.h>
#endif

#include "IMlib.h"

#include "rename.h"
#include "defs.h"
#include "struct.h"
#include "defs.h"
#include "macros.h"
#include "proto.h"
