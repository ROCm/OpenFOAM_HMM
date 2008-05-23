/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * rename.h
 *
 * This file contains renaming #defines to remove any conflicts that the
 * library may have with the users functions.
 *
 * George Irene
 */


/* comm.c */
#define MGridCommInterfaceData MGridCommInterfaceData___
#define MGridCommChangedInterfaceData MGridCommChangedInterfaceData___
#define MGridGlobalSEMax MGridGlobalSEMax___
#define MGridGlobalSEMaxDouble MGridGlobalSEMaxDouble___
#define MGridGlobalSEMin MGridGlobalSEMin___
#define MGridGlobalSESum MGridGlobalSESum___
#define MGridGlobalSESumReal MGridGlobalSESumReal___

/* debug.c */
#define PrintGraph PrintGraph___
#define PrintGraph2 PrintGraph2___

/* grsetup.c */
#define SetUpMGridGraph SetUpMGridGraph___
#define SetUpMGridCtrl SetUpMGridCtrl___
#define ChangeNumbering ChangeNumbering___

/* ikeysort.c */
#define ikeysort ikeysort___

/* memory.c */
#define PreAllocateMGridMemory PreAllocateMGridMemory___
#define FreeMGridWSpace FreeMGridWSpace___
#define FreeMGridCtrl FreeMGridCtrl___
#define CreateMGridGraph CreateMGridGraph___
#define InitMGridGraph InitMGridGraph___
#define FreeMGridGraph FreeMGridGraph___
#define FreeMGridGraphContent FreeMGridGraphContent___

/* move.c */
#define MoveMGridGraph MoveMGridGraph___

/* setup.c */
#define MGridSetUp MGridSetUp___

/* util.c */
#define MGriderrexit MGriderrexit___
#define MGridmyprintf MGridmyprintf___
#define MGridrprintf MGridrprintf___
