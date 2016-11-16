/*---------------------------------------------------------------------------*\
   IH-Cantabria 2015 (http://www.ihcantabria.com/en/)
   IHFOAM 2015 (http://ihfoam.ihcantabria.com/) 

   Author(s):  Javier Lopez Lara (jav.lopez@unican.es)
               Gabriel Barajas   (barajasg@unican.es)
\*---------------------------------------------------------------------------*/

    if ( waveTheory_ == "StokesI" )
    {
            calculatedLevel = RealwaterDepth_ 
	        + timeMult * StokesIFun :: eta 
		(
		    waveHeight_, 
                    waveKx, 
                    xGroup[1],  
                    waveKy, 
                    yGroup[1],  
                    waveOmega, 
                    currTime, 
                    wavePhase_
                );
    }
    else if ( waveTheory_ == "StokesII" )
    {
            calculatedLevel = RealwaterDepth_ 
	        + timeMult * StokesIIFun :: eta 
	        (
                    waveHeight_, 
                    RealwaterDepth_,  
                    waveKx, 
                    xGroup[1],  
                    waveKy, 
                    yGroup[1],  
                    waveOmega, 
                    currTime, 
                    wavePhase_
                );
    }
    else  if ( waveTheory_ == "StokesV" )
    {
            calculatedLevel = RealwaterDepth_ 
	        + timeMult * stokesVFun :: eta 
		(
                    RealwaterDepth_,  
		    waveKx, 
                    waveKy, 
                    lambdaStokesV_, 
                    wavePeriod_, 
                    xGroup[1], 
                    yGroup[1], 
                    currTime, 
                    wavePhase_
                );
    }
    else if ( waveTheory_ == "Cnoidal" )
    {
            calculatedLevel = RealwaterDepth_ 
	        + timeMult * cnoidalFun :: eta 
	        (
		    waveHeight_, 
                    mCnoidal_, 
                    waveKx, 
                    waveKy, 
                    wavePeriod_, 
                    xGroup[1], 
                    yGroup[1], 
                    currTime
                );
    }
    else 
    {
        FatalError << "Wave theory not supported, use:\n" 
	    << "StokesI, StokesII, StokesV, Cnoidal, SolitaryBoussinesq." 
	    << exit(FatalError);
    }
