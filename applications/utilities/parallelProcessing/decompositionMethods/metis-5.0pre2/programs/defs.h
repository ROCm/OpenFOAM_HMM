/*
 * defs.h
 *
 * This file contains various constant definitions
 *
 * Started 8/9/02
 * George
 *
 */

#define CMD_MTYPE               1
#define CMD_ITYPE               2
#define CMD_RTYPE               3

#define CMD_BALANCE             10
#define CMD_NITER               20
#define CMD_NTRIALS             21

#define CMD_TPWGTS              30

#define CMD_SEED                50

#define CMD_OUTPUT              100
#define CMD_NOOUTPUT            101

#define CMD_DBGLVL              1000
#define CMD_HELP                1001




/* The text labels for MTypes */
static char mtypenames[][10] = {"", "None", "MAXTF", "SQRT", "LOG", "IDF"};


/* The text labels for ITypes */
static char itypenames[][10] = {"", "None", "IDF"};


/* The text labels for RTypes */
static char rtypenames[][20]  = {"", "I1", "I2", "E1", "G1", "G1'", "H1", "H2", "SLINK", 
                              "SLINK_W", "CLINK", "CLINK_W", "UPGMA", "UPGMA_W", 
                              "UPGMA_W2", "Cut", "RCut", "NCut", "MMCut"};


