/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * proto.h
 *
 * This file contains header files
 *
 * George Irene
 */

/* coarsen.c */
void LocalMatch_HEM(CtrlType *, GraphType *, WorkSpaceType *);
void Local_CreateCoarseGraph(CtrlType *, GraphType *, WorkSpaceType *, int);


/* comm.c */
void CommInterfaceData(CtrlType *, GraphType *, idxtype *, idxtype *, idxtype *);
void CommChangedInterfaceData(CtrlType *, GraphType *, int, idxtype *, idxtype *, KeyValueType *, KeyValueType *, idxtype *);
int GlobalSEMax(CtrlType *, int);
double GlobalSEMaxDouble(CtrlType *, double);
int GlobalSEMin(CtrlType *, int);
int GlobalSESum(CtrlType *, int);
realtype GlobalSESumReal(CtrlType *, realtype);


/* debug.c */
void PrintVector2(CtrlType *, int, int, idxtype *, char *);


/* diffuse.c */
void ParMETIS_RepartLDiffusion(idxtype *, idxtype *, idxtype *, idxtype *, realtype *, int *, int *, int *, int *, idxtype *, MPI_Comm *);


/* drivers.c */
void AdaptiveUndirected_Partition(CtrlType *, GraphType *, WorkSpaceType *);


/* edge_refine.c */
void ProjectPartition(CtrlType *, GraphType *, WorkSpaceType *);
void ComputePartitionParams(CtrlType *, GraphType *, WorkSpaceType *);
void KWayRefineClean(CtrlType *, GraphType *, WorkSpaceType *, int, float);
void KWayAdaptiveRefineClean(CtrlType *, GraphType *, WorkSpaceType *, int, float);


/* fused.c */
void ParMETIS_FusedElementGraph(idxtype *, idxtype *, realtype *, realtype *, idxtype *, idxtype *,
                       realtype *, int *, int *, int *, int*,  idxtype *, MPI_Comm *);
void CreateFusedElementGraph(CtrlType *, GraphType *, WorkSpaceType *, int *);


/* grsetup.c */
GraphType *SetUpGraph(CtrlType *, idxtype *, idxtype *, idxtype *, idxtype *, realtype *, int);
void SetUpCtrl(CtrlType *ctrl, int, int *, MPI_Comm);
void ChangeNumbering(idxtype *, idxtype *, idxtype *, idxtype *, int, int, int);
void ComputeMoveStatistics(CtrlType *, GraphType *, int *, int *, int *);


/* iidxsort.c */
void iidxsort(int, idxtype *);


/* ikeysort.c */
void ikeysort(int, KeyValueType *);


/* memory.c */
void PreAllocateMemory(CtrlType *, GraphType *, WorkSpaceType *);
void FreeWSpace(WorkSpaceType *);
void FreeCtrl(CtrlType *);
GraphType *CreateGraph(void);
void InitGraph(GraphType *);
void FreeGraph(GraphType *);
void FreeGraphContent(GraphType *);
void FreeInitialGraphAndRemap(GraphType *, int);


/* remap.c */
void ReMapGraph(CtrlType *, GraphType *, int, WorkSpaceType *);
void ComputeTotalVReMap1(CtrlType *, idxtype *, idxtype *, WorkSpaceType *);


/* setup.c */
void SetUp(CtrlType *, GraphType *, WorkSpaceType *);


/* timer.c */
void InitTimers(CtrlType *);
void PrintTimingInfo(CtrlType *);
void PrintTimer(CtrlType *, timer, char *);


/* util.c */
void errexit(char *,...);
void myprintf(CtrlType *, char *f_str,...);
void rprintf(CtrlType *, char *f_str,...);
#ifndef DMALLOC
int *imalloc(int, char *);
idxtype *idxmalloc(int, char *);
float *fmalloc(int, char *);
realtype *realmalloc(int, char *);
int *ismalloc(int, int, char *);
idxtype *idxsmalloc(int, idxtype, char *);
realtype *realsmalloc(int, realtype, char *);
#endif
int *iset(int n, int val, int *x);
idxtype *idxset(int n, idxtype val, idxtype *x);
realtype *realset(int n, realtype val, realtype *x);
int idxamax(int n, idxtype *x);
int idxamin(int n, idxtype *x);
int charsum(int n, char *x);
int isum(int n, int *x);
int idxsum(int, idxtype *);
float snorm2(int, float *);
float sdot(int n, float *, float *);
void saxpy(int, float, float *, float *);
void ikeyvalsort_org(int, KeyValueType *);
int IncKeyValueCmp(const void *, const void *);
void dkeyvalsort(int, KeyValueType *);
int DecKeyValueCmp(const void *, const void *);
int BSearch(int, idxtype *, int);
void RandomPermute(int, idxtype *, int);
void FastRandomPermute(int, idxtype *, int);
int ispow2(int);
int log2(int);



/* var.c */
void ChangeWeights(int, idxtype *, idxtype *, idxtype *, realtype *, MPI_Comm);


/*
void iintsort(int, int *);
void ikeyvalsort(int, KeyValueType *);


double drand48();
void srand48(long);
*/
