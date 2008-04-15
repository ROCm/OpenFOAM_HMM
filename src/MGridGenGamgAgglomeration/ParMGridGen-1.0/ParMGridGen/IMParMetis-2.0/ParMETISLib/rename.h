/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * rename.h
 *
 * This file contains renaming #defines to remove any conflicts that the
 * library may have with the users functions.
 *
 * George
 */

/* coarsen.c */
#define LocalMatch_HEM ___LocalMatch_HEM
#define Local_CreateCoarseGraph ___Local_CreateCoarseGraph

/* comm.c */
#define CommInterfaceData ___CommInterfaceData
#define CommChangedInterfaceData ___CommChangedInterfaceData
#define GlobalSEMax ___GlobalSEMax
#define GlobalSEMaxDouble ___GlobalSEMaxDouble
#define GlobalSEMin ___GlobalSEMin
#define GlobalSESum ___GlobalSESum
#define GlobalSESumReal ___GlobalSESumReal

/* debug.c */
#define PrintVector2 ___PrintVector2

/* diffuse.c */
#define ParMETIS_RepartLDiffusion ___ParMETIS_RepartLDiffusion

/* drivers.c */
#define AdaptiveUndirected_Partition ___AdaptiveUndirected_Partition

/* edge_refine.c */
#define ProjectPartition ___ProjectPartition
#define ComputePartitionParams ___ComputePartitionParams
#define KWayRefineClean ___KWayRefineClean
#define KWayAdaptiveRefineClean ___KWayAdaptiveRefineClean

/* grsetup.c */
#define SetUpGraph ___SetUpGraph
#define SetUpCtrl ___SetUpCtrl
#define ChangeNumbering ___ChangeNumbering
#define ComputeMoveStatistics ___ComputeMoveStatistics

/* iidxsort.c */
#define iidxsort ___iidxsort

/* ikeysort.c */
#define ikeysort ___ikeysort

/* memory.c */
#define PreAllocateMemory ___PreAllocateMemory
#define FreeWSpace ___FreeWSpace
#define FreeCtrl ___FreeCtrl
#define CreateGraph ___CreateGraph
#define InitGraph ___InitGraph
#define FreeGraph ___FreeGraph
#define FreeGraphContent ___FreeGraphContent
#define FreeInitialGraphAndRemap ___FreeInitialGraphAndRemap

/* remap.c */
#define ReMapGraph ___ReMapGraph
#define ComputeTotalVReMap1 ___ComputeTotalVReMap1

/* setup.c */
#define SetUp ___SetUp

/* timer.c */
#define InitTimers ___InitTimers
#define PrintTimingInfo ___PrintTimingInfo
#define PrintTimer ___PrintTimer

/* util.c */
#define errexit ___errexit
#define myprintf ___myprintf
#define rprintf ___rprintf
#ifndef DMALLOC
#define imalloc ___imalloc
#define idxmalloc ___idxmalloc
#define fmalloc ___fmalloc
#define realmalloc ___realmalloc
#define ismalloc ___ismalloc
#define idxsmalloc ___idxsmalloc
#define realsmalloc ___realsmalloc
#endif
#define iset ___iset
#define idxset ___idxset
#define realset ___realset
#define idxamax ___idxamax
#define idxamin ___idxamin
#define charsum ___charsum
#define isum ___isum
#define idxsum ___idxsum
#define snorm2 ___snorm2
#define sdot ___sdot
#define saxpy ___saxpy
#define ikeyvalsort_org ___ikeyvalsort_org
#define IncKeyValueCmp ___IncKeyValueCmp
#define dkeyvalsort ___dkeyvalsort
#define DecKeyValueCmp ___DecKeyValueCmp
#define BSearch ___BSearch
#define RandomPermute ___RandomPermute
#define FastRandomPermute ___FastRandomPermute
#define ispow2 ___ispow2
#define log2 ___log2
