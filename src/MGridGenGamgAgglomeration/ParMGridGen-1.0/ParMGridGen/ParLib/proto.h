/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * proto.h
 *
 * This file contains header files
 *
 * George Irene
 */

/* parmgridgen.c */
void ParMGridGen(idxtype *, idxtype *, realtype *, realtype *,
      idxtype *, realtype *, int *, int, int, int *, idxtype *,
      MPI_Comm *);
void MoveRefine(idxtype *, idxtype **, realtype **, realtype **,
      idxtype **, realtype **, idxtype **, int *, int, int, int *,
      int *, int *, idxtype **, MPI_Comm *);
void PrintAspectRatioStats(idxtype *, idxtype *, realtype *, realtype *,
      idxtype *, realtype *, int *, int, int, int *, idxtype *, MPI_Comm);


/* comm.c */
void MGridCommInterfaceData(MGridCtrlType *, MGridGraphType *, idxtype *,
      idxtype *, idxtype *);
void MGridCommChangedInterfaceData(MGridCtrlType *, MGridGraphType *, int, 
      idxtype *, idxtype *, KeyValueType *, KeyValueType *, idxtype *);
int MGridGlobalSEMax(MGridCtrlType *, int);
double MGridGlobalSEMaxDouble(MGridCtrlType *, double);
int MGridGlobalSEMin(MGridCtrlType *, int);
int MGridGlobalSESum(MGridCtrlType *, int);
realtype MGridGlobalSESumReal(MGridCtrlType *, realtype);


/* debug.c */
void PrintGraph(MGridCtrlType *, MGridGraphType *);
void PrintGraph2(MGridCtrlType *, MGridGraphType *);
void WriteMetisGraph(int, idxtype *, idxtype *, idxtype *, realtype *);
void WriteMetisGraph2(int, idxtype *, idxtype *, idxtype *, realtype *, MPI_Comm);


/* grsetup.c */
MGridGraphType *SetUpMGridGraph(MGridCtrlType *, idxtype *, idxtype *,
      realtype *, realtype *, idxtype *, realtype *);
void SetUpMGridCtrl(MGridCtrlType *ctrl, int, int, int *, MPI_Comm);
void ChangeNumbering(idxtype *, idxtype *, idxtype *, idxtype *, int, int, int);
void RemoveEdges(idxtype *, idxtype *, idxtype *, realtype *,  idxtype *,
      idxtype *, realtype *, MPI_Comm);
void CorrectSurfaceEdges(idxtype *, idxtype *, idxtype *, realtype *, 
      realtype *, MPI_Comm);
void RemoveCorrectSurfaceEdges(idxtype *, idxtype *, idxtype *, realtype *,
      idxtype *, idxtype *, realtype *, realtype *, MPI_Comm);


/* ikeysort.c */
void ikeysort(int, KeyValueType *);


/* memory.c */
void PreAllocateMGridMemory(MGridCtrlType *, MGridGraphType *, MGridWorkSpaceType *);
void FreeMGridWSpace(MGridWorkSpaceType *);
void FreeMGridCtrl(MGridCtrlType *);
MGridGraphType *CreateMGridGraph(void);
void InitMGridGraph(MGridGraphType *);
void FreeMGridGraph(MGridGraphType *);
void FreeMGridGraphContent(MGridGraphType *);


/* move.c */
MGridGraphType *MoveMGridGraph(MGridCtrlType *, MGridGraphType *,
      MGridWorkSpaceType *);


/* setup.c */
void MGridSetUp(MGridCtrlType *, MGridGraphType *, MGridWorkSpaceType *);


/* util.c */
void MGriderrexit(char *,...);
void MGridmyprintf(MGridCtrlType *, char *f_str,...);
void MGridrprintf(MGridCtrlType *, char *f_str,...);


/* parmgridgen.c */
void TestParMGridGen(char *, int *, int, int, MPI_Comm);


/* tio.c */
void MGridReadTestGraph(MGridGraphType *, char *, MPI_Comm);
double *ReadTestCoordinates(MGridGraphType *, char *, int, MPI_Comm);
void ReadMGridGraph(char *, int *, idxtype **, idxtype **, realtype **, realtype **, realtype **);
void WriteParallelPartition(char *, idxtype *, idxtype *, int, int, int);  

