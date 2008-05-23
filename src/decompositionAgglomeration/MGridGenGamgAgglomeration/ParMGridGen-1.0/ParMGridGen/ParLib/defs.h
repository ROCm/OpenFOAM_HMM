/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * defs.h
 *
 * This file contains constant definitions
 *
 * George Irene
 */

#define MAXLINE                  1280000

/* CType Schemes */
#define MATCH_RM                1
#define MATCH_HEM               2
#define MATCH_HEM_SLOW          3
#define MATCH_HEM_TRUE          4

/* RType Schemes */
#define REFINE_AR               1
#define REFINE_WAR              2
#define REFINE_SCUT             3
#define REFINE_MINMAXAVAR       4
#define REFINE_MINMAXAR         5
#define REFINE_MULTIOBJECTIVE   6
#define REFINE_MULTIOBJECTIVE2  7

#define NGD_PASSES		20

/* Meaning of various options[] parameters 3-6 */
#define OPTION_DBGLVL		2
#define OPTION_CTYPE            0
#define OPTION_RTYPE            1
#define OPTION_DIM            	3

#define XYZ_XCOORD		1
#define XYZ_SPFILL		2

/* Type of initial vertex separator algorithms */
#define ISEP_EDGE		1
#define ISEP_NODE		2

#define UNMATCHED		-1
#define MAYBE_MATCHED		-2
#define TOO_HEAVY		-3

#define MAX_PES			8192
#define MAX_INT			(1<<30)

#define HTABLE_EMPTY    	-1

#define NGR_PASSES		4	/* Number of greedy refinement passes */
#define NLGR_PASSES		5	/* Number of GR refinement during IPartition */

#define IGNORE_FRACTION		0.9	/* What fraction of vertices will not be colored */

#define KEEP_BIT        (idxtype)536870912        /* 1<<29 */

#define UNBALANCE_FRACTION		1.05
#define ORDER_UNBALANCE_FRACTION	1.05

#define MAXVWGT_FACTOR		1.4

#define MATCH_LOCAL		1
#define MATCH_GLOBAL		2

#define EDGE_WEIGHT		1
#define NOEDGE_WEIGHT		2

/* Debug Levels */
#define DBG_TRACK	128 		/* Track flow of program Irene Debugging */
#define DBG_STEPS	256 		/* Track convergence through iterations */

