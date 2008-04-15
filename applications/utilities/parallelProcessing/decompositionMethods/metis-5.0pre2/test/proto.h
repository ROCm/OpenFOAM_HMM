/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * proto.h
 *
 * This file contains header files
 *
 * Started 10/19/95
 * George
 *
 * $Id: proto.h 1428 2007-04-06 23:37:27Z karypis $
 *
 */

#ifndef _TEST_PROTO_H_
#define _TEST_PROTO_H_

void ReadGraph(GraphType *, char *, idxtype *);

void Test_PartGraph(idxtype, idxtype *, idxtype *);
int  VerifyPart(idxtype, idxtype *, idxtype *, idxtype *, idxtype *, idxtype, idxtype, idxtype *);
int  VerifyWPart(idxtype, idxtype *, idxtype *, idxtype *, idxtype *, idxtype, float *, idxtype, idxtype *);
void Test_PartGraphV(idxtype, idxtype *, idxtype *);
int  VerifyPartV(idxtype, idxtype *, idxtype *, idxtype *, idxtype *, idxtype, idxtype, idxtype *);
int  VerifyWPartV(idxtype, idxtype *, idxtype *, idxtype *, idxtype *, idxtype, float *, idxtype, idxtype *);
void Test_PartGraphmC(idxtype, idxtype *, idxtype *);
int  VerifyPartmC(idxtype, idxtype, idxtype *, idxtype *, idxtype *, idxtype *, idxtype, float *, idxtype, idxtype *);
void Test_ND(idxtype, idxtype *, idxtype *);
int  VerifyND(idxtype, idxtype *, idxtype *);

#endif
