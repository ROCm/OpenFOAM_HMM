/*---------------------------------------------------------------------------*\
   IH-Cantabria 2015 (http://www.ihcantabria.com/en/)
   IHFOAM 2015 (http://ihfoam.ihcantabria.com/) 

   Author(s):  Javier Lopez Lara (jav.lopez@unican.es)
               Gabriel Barajas   (barajasg@unican.es)
\*---------------------------------------------------------------------------*/

    if ( waveTheory_ == "StokesI" )
    {
        forAll(calculatedLevel, itS1)
        {
            calculatedLevel[itS1] = RealwaterDepth_ 
		+ timeMult * StokesIFun :: eta 
		(
    		    waveHeight_, 
	            waveKx, 
	            xGroup[itS1], 
                    waveKy, 
                    yGroup[itS1], 
                    waveOmega, 
                    currTime, 
                    wavePhase_
                );
        }
    }
    else if ( waveTheory_ == "StokesII" )
    {
        forAll(calculatedLevel, itS2)
        {
            calculatedLevel[itS2] = RealwaterDepth_ 
	        + timeMult * StokesIIFun :: eta 
		(
                    waveHeight_, 
                    RealwaterDepth_,  
                    waveKx, 
                    xGroup[itS2], 
                    waveKy, 
                    yGroup[itS2], 
                    waveOmega, 
                    currTime, 
                    wavePhase_
                );
        }
    }
    else  if ( waveTheory_ == "StokesV" )
    {
        forAll(calculatedLevel, it1)
        {
            calculatedLevel[it1] = RealwaterDepth_ 
	        + timeMult * stokesVFun :: eta 
		(
                    RealwaterDepth_,  
		    waveKx, 
                    waveKy, 
                    lambdaStokesV_, 
                    wavePeriod_, 
                    xGroup[it1], 
                    yGroup[it1], 
                    currTime, 
                    wavePhase_
                );
        }
    }
    else if ( waveTheory_ == "Cnoidal" )
    {
        forAll(calculatedLevel, it3)
        {
            calculatedLevel[it3] = RealwaterDepth_ 
	        + timeMult * cnoidalFun :: eta 
	        (
		    waveHeight_, 
                    mCnoidal_, 
                    waveKx, 
                    waveKy, 
                    wavePeriod_, 
                    xGroup[it3], 
                    yGroup[it3], 
                    currTime
                );
        }
    }
    else 
    {
        FatalError << "Wave theory not supported, use:\n" 
	    << "StokesI, StokesII, StokesV, Cnoidal, SolitaryBoussinesq" 
	    << exit(FatalError);
    }
