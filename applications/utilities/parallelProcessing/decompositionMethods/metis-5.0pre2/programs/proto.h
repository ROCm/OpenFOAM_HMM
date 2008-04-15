/*
 * proto.h 
 *
 * This file contains function prototypes
 *
 * Started 11/1/99
 * George
 *
 * $Id: proto.h,v 1.15 2003/04/04 23:22:49 karypis Exp $
 *
 */

#ifndef _PROTOBIN_H_
#define _PROTOBIN_H_


/* io.c */ 
void ReadGraph(GraphType *, char *, idxtype *); 
void ReadCoordinates(GraphType *, char *); 
void WritePartition(char *, idxtype *, idxtype, idxtype);
void WriteMeshPartition(char *, idxtype, idxtype, idxtype *, idxtype, idxtype *);
void WritePermutation(char *, idxtype *, idxtype);
int CheckGraph(GraphType *);
int MeshType(char *);
idxtype *ReadWgt(char *, idxtype *, idxtype *, idxtype *); 
idxtype *ReadMesh(char *, idxtype *, idxtype *, idxtype *); 
idxtype *ReadMeshWgt(char *, idxtype *, idxtype *, idxtype *, idxtype *);
idxtype *ReadMixedMesh(char *, idxtype *, idxtype *, idxtype *);
idxtype *ReadMixedMeshWgt(char *, idxtype *, idxtype *, idxtype *, idxtype *);
void WriteGraph(char *, idxtype, idxtype *, idxtype *);
idxtype MixedElements(char *);
idxtype *ReadMgcnums(char *);
void WriteWgtGraph(char *, idxtype , idxtype *, idxtype *, idxtype *);


/* smbfactor.c */
void ComputeFillIn(GraphType *, idxtype *);
idxtype ComputeFillIn2(GraphType *, idxtype *);
idxtype smbfct(idxtype, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *); 

/* cmdline.c */
void parse_cmdline(ParamType *params, int argc, char *argv[]);


/* cpmetis.c */

#endif 
