/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * aratio.c
 *
 * This file contains the aspect ratio definition for 2D and 3D
 *
 * Irene
 */

#include "mgridgen.h"

ASPECTRATIOFUNCTION ARATIO1;
ASPECTRATIOFUNCTION ARATIO;
ASPECTRATIOFUNCTION ARATIO2;


/*************************************************************************
* Aspect Ratio Definition for 2-D problems         
*************************************************************************/
realtype ARATIO1_2D(realtype circumf, realtype surf)
{
realtype ar;

ar = pow((circumf), 2)/(surf);

return ar;
}


realtype ARATIO_2D(realtype circumf, realtype surf)
{
realtype ar;

ar = ((circumf)*(circumf))/(surf);

return ar;
}


realtype ARATIO2_2D(realtype circumf, realtype surf)
{
realtype ar;

ar = ((circumf)*(circumf)*(circumf)*(circumf))/((surf)*(surf));

return ar;
}


/*************************************************************************
* Aspect Ratio Definition for 3-D problems         
*************************************************************************/
realtype ARATIO1_3D(realtype surf, realtype vol)
{
realtype ar;

ar = pow((surf), 1.5)/(vol);

return ar;
}


realtype ARATIO_3D(realtype surf, realtype vol)
{
realtype ar;

ar = sqrt((surf)*(surf)*(surf))/(vol);

return ar;
}


realtype ARATIO2_3D(realtype surf, realtype vol)
{
realtype ar;

ar = ((surf)*(surf)*(surf))/((vol)*(vol));

return ar;
}
