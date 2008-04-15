/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * defs.h
 *
 * This file contains constant definitions
 *
 * George Irene
 */

#define LTERM                   (void **) 0     /* List terminator for GKfree() */

#define NGD_PASSES		20

#define OPTION_IPART		1
#define OPTION_FOLDF		2
#define OPTION_DBGLVL		3

#define XYZ_XCOORD		1
#define XYZ_SPFILL		2

/* Type of initial partitioning algorithms */
#define IPART_SER		1
#define IPART_RB		2

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
#define NIPARTS			8	/* Number of random initial partitions */
#define NLGR_PASSES		5	/* Number of GR refinement during IPartition */

#define IGNORE_FRACTION		0.9	/* What fraction of vertices will not be colored */

#define KEEP_BIT        (idxtype)536870912        /* 1<<29 */

#define COARSEN_FRACTION	0.75	/* Node reduction between succesive coarsening levels */
#define COARSEN_FRACTION2	0.55	/* Node reduction between succesive coarsening levels */
#define UNBALANCE_FRACTION		1.05
#define ORDER_UNBALANCE_FRACTION	1.05

#define MAXVWGT_FACTOR		1.4

#define MATCH_LOCAL		1
#define MATCH_GLOBAL		2

#define EDGE_WEIGHT		1
#define NOEDGE_WEIGHT		2

/* Debug Levels */
#define DBG_TIME	1		/* Perform timing analysis */
#define DBG_INFO	2		/* Perform timing analysis */
#define DBG_PROGRESS   	4		/* Show the coarsening progress */
#define DBG_REFINEINFO	8		/* Show info on communication during folding */
#define DBG_MATCHINFO	16		/* Show info on matching */
#define DBG_RMOVEINFO	32		/* Show info on communication during folding */
#define DBG_REMAP	64		/* Determines if remapping will take place */
#define DBG_TRACK	128		/* IRENE */
