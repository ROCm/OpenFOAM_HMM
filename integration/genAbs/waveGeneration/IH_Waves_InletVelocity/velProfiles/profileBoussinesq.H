/*---------------------------------------------------------------------------*\
   IH-Cantabria 2015 (http://www.ihcantabria.com/en/)
   IHFOAM 2015 (http://ihfoam.ihcantabria.com/) 

   Author(s):  Javier Lopez Lara (jav.lopez@unican.es)
               Gabriel Barajas   (barajasg@unican.es)
\*---------------------------------------------------------------------------*/

    if (zSup[cellIndex] < heights[cellGroup[cellIndex]-1])
    {
        patchU[cellIndex] = BoussinesqFun :: U 
	    (
                waveHeight_, 
                RealwaterDepth_, 
                xGroup[cellGroup[cellIndex]-1], 
                yGroup[cellGroup[cellIndex]-1], 
                waveAngle, 
                currTime, 
                X0, 
                patchHeight[cellIndex]
             );

        patchV[cellIndex] = patchU[cellIndex] * sin(waveAngle) * timeMult;
        patchU[cellIndex] = patchU[cellIndex] * cos(waveAngle) * timeMult;

        patchW[cellIndex] = BoussinesqFun :: W 
	    (
                waveHeight_, 
                RealwaterDepth_, 
                xGroup[cellGroup[cellIndex]-1], 
                yGroup[cellGroup[cellIndex]-1], 
                waveAngle, 
                currTime, 
                X0, 
                patchHeight[cellIndex]
            );

	patchW[cellIndex] = patchW[cellIndex] * timeMult;

    }
    else if ( zInf[cellIndex] > calculatedLevel[cellGroup[cellIndex]-1])
    {
        patchU[cellIndex] = 0.0;
        patchV[cellIndex] = 0.0;
        patchW[cellIndex] = 0.0;
    } 
    else 
    {	
        if ( calculatedLevel[cellGroup[cellIndex]-1] 
	    < measuredLevels[cellGroup[cellIndex]-1] )
        {
	    if ((zSup[cellIndex] > calculatedLevel[cellGroup[cellIndex]-1]) 
	        & (zInf[cellIndex] < calculatedLevel[cellGroup[cellIndex]-1]))
	    {
                auxiliar = calculatedLevel[cellGroup[cellIndex]-1] 
		    - zInf[cellIndex];
		auxiliarTotal = zSup[cellIndex]-zInf[cellIndex];
		auxiliarTotal = auxiliar/auxiliarTotal;

		patchU[cellIndex] = BoussinesqFun :: U 
		    (
                        waveHeight_, 
                    	RealwaterDepth_, 
                    	xGroup[cellGroup[cellIndex]-1], 
                    	yGroup[cellGroup[cellIndex]-1], 
                    	waveAngle, 
                    	currTime, 
                    	X0, 
                    	patchHeight[cellIndex]
                    );

		patchV[cellIndex] = auxiliarTotal * patchU[cellIndex] 
		    * sin(waveAngle) * timeMult;
		patchU[cellIndex] = auxiliarTotal * patchU[cellIndex] 
		    * cos(waveAngle) * timeMult;

		patchW[cellIndex] = BoussinesqFun :: W 
		    (
                        waveHeight_, 
                    	RealwaterDepth_, 
                    	xGroup[cellGroup[cellIndex]-1], 
                    	yGroup[cellGroup[cellIndex]-1], 
                    	waveAngle, 
                    	currTime, 
                    	X0, 
                    	patchHeight[cellIndex]
                    );

		patchW[cellIndex] = auxiliarTotal * patchW[cellIndex]
		    * timeMult;	
	    }
	}	
	else if ( calculatedLevel[cellGroup[cellIndex]-1] 
	    > measuredLevels[cellGroup[cellIndex]-1] )
	{
	    if (zSup[cellIndex] < calculatedLevel[cellGroup[cellIndex]-1])
	    {
	        patchU[cellIndex] = BoussinesqFun :: U 
		    (
                        waveHeight_, 
                        RealwaterDepth_, 
                        xGroup[cellGroup[cellIndex]-1], 
                        yGroup[cellGroup[cellIndex]-1], 
                        waveAngle, 
                        currTime, 
                        X0, 
                        measuredLevels[cellGroup[cellIndex]-1]
                    );

		patchV[cellIndex] = patchU[cellIndex] * sin(waveAngle) 
		    * timeMult;
		patchU[cellIndex] = patchU[cellIndex] * cos(waveAngle) 
		    * timeMult;

		patchW[cellIndex] = BoussinesqFun :: W 
		    (
                        waveHeight_, 
                    	RealwaterDepth_, 
                    	xGroup[cellGroup[cellIndex]-1], 
                    	yGroup[cellGroup[cellIndex]-1], 
                    	waveAngle, 
                    	currTime, 
                    	X0, 
                    	measuredLevels[cellGroup[cellIndex]-1]
                    );

		patchW[cellIndex] = patchW[cellIndex]*timeMult;			
	    }
	    else if ((zSup[cellIndex] > calculatedLevel[cellGroup[cellIndex]-1]) 
	        & (zInf[cellIndex] < calculatedLevel[cellGroup[cellIndex]-1]) )
	    {
	        auxiliar = calculatedLevel[cellGroup[cellIndex]-1] 
		     - zInf[cellIndex];
		auxiliarTotal = zSup[cellIndex]-zInf[cellIndex];
		auxiliarTotal = auxiliar/auxiliarTotal;

		patchU[cellIndex] = BoussinesqFun :: U 
		    (
		        waveHeight_, 
			RealwaterDepth_, 
			xGroup[cellGroup[cellIndex]-1], 
			yGroup[cellGroup[cellIndex]-1], 
			waveAngle, 
			currTime, 
			X0, 
			measuredLevels[cellGroup[cellIndex]-1]
		    );

		patchV[cellIndex] = auxiliarTotal * patchU[cellIndex] 
		    * sin(waveAngle) * timeMult;
		patchU[cellIndex] = auxiliarTotal * patchU[cellIndex] 
		    * cos(waveAngle) * timeMult;

		patchW[cellIndex] = BoussinesqFun :: W 
		    (
                        waveHeight_, 
                    	RealwaterDepth_, 
                    	xGroup[cellGroup[cellIndex]-1], 
                    	yGroup[cellGroup[cellIndex]-1], 
                    	waveAngle, 
                    	currTime, 
                    	X0, 
                    	measuredLevels[cellGroup[cellIndex]-1]
                   );

		patchW[cellIndex] = auxiliarTotal * patchW[cellIndex] 
		    * timeMult;	
            }
        }	
    }
