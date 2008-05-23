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
GraphType *Coarsen(CtrlType *, GraphType *);
GraphType *Coarsen_Restricted(CtrlType *, GraphType *);
void CreateCoarseGraph(GraphType *, int, idxtype *, idxtype *);
GraphType *SetUpCoarseGraph(GraphType *, int);


/* io.c */
void ReadGraph(GraphType *, char *);
void TransformGraph(char *);
void WritePartition(char *, idxtype *, int, int);
void PrintGraph(GraphType *);


/* kwayfm.c */
void Random_KWayARatioRefine(CtrlType *, GraphType *, int);
void Random_KWayWeightARatioRefine(CtrlType *, GraphType *, int);
void Random_KWaySCutRefine(CtrlType *, GraphType *, int);
void Random_KWayMinMaxAverageARatioRefine(CtrlType *, GraphType *, int);
void Random_KWayMinMaxARatioRefine(CtrlType *, GraphType *, int);
void Random_KWayMultiObjRefine(CtrlType *, GraphType *, int);
void Random_KWayMultiObjRefine2(CtrlType *, GraphType *, int);
void CheckParams(CtrlType *, GraphType *);


/* match.c */
void Match_RM(CtrlType *, GraphType *);
void Match_HEM(CtrlType *, GraphType *);
void Match_HEM_Slow(CtrlType *, GraphType *);
void Match_HEM_Slow_Restricted(CtrlType *ctrl, GraphType *graph);
void Match_HEM_True(CtrlType *, GraphType *);


/* merge.c */
void Cycle(CtrlType *, GraphType *, int);
void Merge(CtrlType *, GraphType *, int);
void Contribute(CtrlType *, GraphType *, int);
void Merge_ARatio(CtrlType *, GraphType *);
void Merge_WeightARatio(CtrlType *, GraphType *);
void Merge_MinMaxARatio(CtrlType *, GraphType *);
void Merge_MultiObj(CtrlType *, GraphType *);
void Merge0(CtrlType *, GraphType *);
void Contribute_ARatio(CtrlType *, GraphType *);
void Contribute_WeightARatio(CtrlType *, GraphType *);
void Contribute_MinMaxARatio(CtrlType *, GraphType *);
void Contribute_MultiObj(CtrlType *, GraphType *);


/* mgridgen.c */
void MGridGen(int, idxtype *, realtype *, realtype *, idxtype *, realtype *,
              int, int, int *, int *, int *, idxtype *);
void MGridGenRefine(int, idxtype *, realtype *, realtype *, idxtype *,
                    idxtype *, realtype *, int , int , int *, int *, int *,
                    idxtype *);
void CreateGrid(CtrlType *, GraphType *);


/* refine.c */
void RefineKWay(CtrlType *, GraphType *, GraphType *, int);
void RefineKWayOnce(CtrlType *, GraphType *, int);
void ComputeKWayPartitionParams(CtrlType *, GraphType *);
void ProjectKWayPartition(GraphType *);
void ComputeGridStatistics(CtrlType *, GraphType *);
void BreakComponents(CtrlType *, GraphType *);
realtype ComputeFunction(int, CtrlType *, GraphType *);
void ComputeAllFunctions(CtrlType *, GraphType *);


/* setup.c */
void SetUpGraph(GraphType *, int, idxtype *, realtype *, realtype *, idxtype *,
                realtype *);
void FreeGraph(GraphType *);
GraphType *CreateGraph(void);


/* aratio.c */
realtype ARATIO1_2D(realtype, realtype);
realtype ARATIO_2D(realtype, realtype);
realtype ARATIO2_2D(realtype, realtype);
realtype ARATIO1_3D(realtype, realtype);
realtype ARATIO_3D(realtype, realtype);
realtype ARATIO2_3D(realtype, realtype);
typedef realtype (*ASPECTRATIOFUNCTION) (realtype, realtype);
extern ASPECTRATIOFUNCTION ARATIO1;
extern ASPECTRATIOFUNCTION ARATIO;
extern ASPECTRATIOFUNCTION ARATIO2;

