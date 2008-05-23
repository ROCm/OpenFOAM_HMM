/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * defs.h
 *
 * This file contains constant definitions
 *
 * George Irene
 */

#define MAXLINE			1280000


/* Meaning of various options[] parameters */
#define OPTION_CTYPE		0
#define OPTION_RTYPE		1
#define OPTION_DBGLVL		2
#define OPTION_DIM		3


/* CType Schemes */
#define MATCH_RM		1
#define MATCH_HEM		2
#define MATCH_HEM_SLOW          3
#define MATCH_HEM_TRUE          4

/* RType Schemes */
#define REFINE_AR		1
#define REFINE_WAR		2
#define REFINE_SCUT             3
#define REFINE_MINMAXAVAR       4
#define REFINE_MINMAXAR         5
#define REFINE_MULTIOBJECTIVE   6
#define REFINE_MULTIOBJECTIVE2  7


#define UNMATCHED		-1

#define MAXIDX			(1<<30)

#define HTABLE_EMPTY    	-1


/* Debug Levels */
#define DBG_PROGRESS   	1		/* If you should show size info per level */
#define DBG_OUTPUT	2		/* If you should output intermediate partitions */
#define DBG_COARSEN	4
#define DBG_REFINE	8
#define DBG_MOVEINFO	16
#define DBG_MERGE	32		/* If you should check merging */
#define DBG_CONTR	64		/* If you should check contributing */
#define DBG_TRACK	128		/* Track progress; Print Final Results only */
