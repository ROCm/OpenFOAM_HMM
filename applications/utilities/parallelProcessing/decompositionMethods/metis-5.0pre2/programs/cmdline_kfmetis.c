/*
 * cmdline_kfmetis.c
 *
 * This file parses the command line arguments
 *
 * Started 8/9/02
 * George
 *
 * $Id: cmdline_kfmetis.c,v 1.1 2002/08/12 15:20:50 karypis Exp $
 *
 */


#include <metisbin.h>


/*-------------------------------------------------------------------
 * Command-line options 
 *-------------------------------------------------------------------*/
static struct gk_option long_options[] = {
  {"mtype",          1,      0,      CMD_MTYPE},
  {"itype",          1,      0,      CMD_ITYPE},
  {"rtype",          1,      0,      CMD_RTYPE},

  {"balanced",       0,      0,      CMD_BALANCE},

  {"niter",          1,      0,      CMD_NITER},

  {"tpwgts",         1,      0,      CMD_TPWGTS},

  {"seed",           1,      0,      CMD_SEED},

  {"dbglvl",         1,      0,      CMD_DBGLVL},

  {"help",           0,      0,      CMD_HELP},
  {0,                0,      0,      0}
};



/*-------------------------------------------------------------------
 * Mappings for the various parameter values
 *-------------------------------------------------------------------*/
static gk_StringMap_t mtype_options[] = {
 {"rm",                 MTYPE_RM},
 {"hem",                MTYPE_HEM},
 {"shem",               MTYPE_SHEM},
 {"shebm",              MTYPE_SHEBM_ONENORM},
 {"sbhem",              MTYPE_SBHEM_ONENORM},
 {NULL,                 0}
};


static gk_StringMap_t itype_options[] = {
 {"greedy",              ITYPE_GGPKL},
 {"random",              ITYPE_RANDOM},
 {NULL,                 0}
};


static gk_StringMap_t rtype_options[] = {
 {"fm",            RTYPE_FM},
 {NULL,                 0}
};



/*-------------------------------------------------------------------
 * Mini help
 *-------------------------------------------------------------------*/
static char helpstr[][100] =
{
" ",
"Usage: pmetis [options] <filename> <nparts>",
" ",
" Required parameters",
"    filename    Stores the graph to be partitioned.",
"    nparts      The number of partitions to split the graph.",
" ",
" Optional parameters",
"  -mtyep=string",
"     Specifies the scheme to be used to match the vertices of the graph",
"     during the coarsening.",
"     The possible values are:",
"        rm       - Random matching",
"        hem      - Heavy-edge matching",
"        shem     - Sorted heavy-edge matching [default]",
"        shebm    - Combination of shem and balanced matching for",
"                   multi-constraint.",
"        sbhem    - Similar as shebm but priority is given to balance",
" ",
"  -itype=string",
"     Specifies the scheme to be used to compute the initial partitioning",
"     of the graph.",
"     The possible values are:",
"        greedy   - Grow a bisection using a greedy strategy [default]",
"        random   - Compute a bisection at random",
" ",
"  -rtype=string",
"     Specifies the scheme to be used for refinement",
"     The possible values are:",
"        fm       - FM refinement",
" ",
"  -balance",
"     Specifies that the final partitioning should contain nparts-1 equal",
"     size partitions with the last partition having upto nparts-1 fewer",
"     vertices.",
" ",
"  -seed=int      ",
"     Selects the seed of the random number generator.  ",
" ",
"  -dbglvl=int      ",
"     Selects the dbglvl.  ",
" ",
"  -help",
"     Prints this message.",
""
};

static char shorthelpstr[][100] = {
" ",
"   Usage: pmetis [options] <filename> <nparts>",
"          use 'pmetis -help' for a summary of the options.",
""
};
 


/*************************************************************************
* This is the entry point of the command-line argument parser
**************************************************************************/
void parse_cmdline(ParamType *params, int argc, char *argv[])
{
  int i, j, k;
  int c, option_index;

  /* initialize the params data structure */
  params->mtype         = PMETIS_CTYPE;
  params->itype         = PMETIS_ITYPE;
  params->rtype         = PMETIS_RTYPE;
  params->dbglvl        = PMETIS_DBGLVL;

  params->balance       = 0;
  params->seed            = -1;
  params->dbglvl          = 0;
  params->filename        = NULL;
  params->nparts          = 1;


  /* Parse the command line arguments  */
  while ((c = gk_getopt_long_only(argc, argv, "", long_options, &option_index)) != -1) {
    switch (c) {
      case CMD_MTYPE:
        if (gk_optarg)
          if ((params->mtype = gk_GetStringID(mtype_options, gk_optarg)) == -1)
            errexit("Invalid option -%s=%s\n", long_options[option_index].name, gk_optarg);
        break;
      case CMD_ITYPE:
        if (gk_optarg)
          if ((params->itype = gk_GetStringID(itype_options, gk_optarg)) == -1)
            errexit("Invalid option -%s=%s\n", long_options[option_index].name, gk_optarg);
        break;
      case CMD_RTYPE:
        if (gk_optarg)
          if ((params->rtype = gk_GetStringID(rtype_options, gk_optarg)) == -1)
            errexit("Invalid option -%s=%s\n", long_options[option_index].name, gk_optarg);
        break;

      case CMD_BALANCE:
        params->balance = 1;
        break;


      case CMD_SEED:
        if (gk_optarg) params->seed = atoi(gk_optarg);
        break;

      case CMD_DBGLVL:
        if (gk_optarg) params->dbglvl = atoi(gk_optarg);
        break;

      case CMD_HELP:
        for (i=0; strlen(helpstr[i]) > 0; i++)
          mprintf("%s\n", helpstr[i]);
        exit(0);
        break;
      case '?':
      default:
        mprintf("Illegal command-line option(s)\nUse %s -help for a summary of the options.\n", argv[0]);
        exit(0);
    }
  }

  if (argc-gk_optind != 2) {
    mprintf("Missing parameters.");
    for (i=0; strlen(shorthelpstr[i]) > 0; i++)
      mprintf("%s\n", shorthelpstr[i]);
    exit(0);
  }

  params->filename = strdup(argv[gk_optind++]);
  params->nparts   = atoi(argv[gk_optind++]);
    
}


