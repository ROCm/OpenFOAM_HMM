/*---------------------------------------------------------------------------*\
   IH-Cantabria 2015 (http://www.ihcantabria.com/en/)
   IHFOAM 2015 (http://ihfoam.ihcantabria.com/) 

   Author(s):  Javier Lopez Lara (jav.lopez@unican.es)
               Gabriel Barajas   (barajasg@unican.es)
\*---------------------------------------------------------------------------*/

    if ( genAbs_ )
    {  
	if ( zInf[cellIndex] >= measuredLevelsGENAB[cellGroup[cellIndex]-1] )
	{
	    patchUABS[cellIndex] = 0.0;
	}
	else 
	{
	    patchUABS[cellIndex] = ( calculatedLevel[cellGroup[cellIndex]-1] - measuredLevelsGENAB[cellGroup[cellIndex]-1] ) * sqrt(g/measuredLevelsGENAB[cellGroup[cellIndex]-1]);
        }
    }
    patchU[cellIndex] = patchU[cellIndex] + leftORright_ * patchUABS[cellIndex];
