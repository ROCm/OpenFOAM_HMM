/*---------------------------------------------------------------------------*\
   IH-Cantabria 2015 (http://www.ihcantabria.com/en/)
   IHFOAM 2015 (http://ihfoam.ihcantabria.com/) 

   Author(s):  Javier Lopez Lara (jav.lopez@unican.es)
               Gabriel Barajas   (barajasg@unican.es)
\*---------------------------------------------------------------------------*/

    if (nPaddles_ < 1)
    {
        FatalError << "Check nPaddles value." << exit(FatalError);
    }

    if ( nPaddles_ > 1 )
    {
        nPaddles_ = decreaseNPaddles( nPaddles_, patchD, dMin, dSpan );
        reduce(nPaddles_, minOp<label>());
    }

    if ( waveHeight_ < 0.0 )
    {
        FatalError << "Check wave height value." << exit(FatalError);
    }

    if ( waveTheory_ == "aaa" )
    {
        FatalError << "Wave theory not specified." << exit(FatalError);
    }
    else if ( waveTheory_ == "Boussinesq" )
    {
        Info << "\nIH Wave Generation BC" << endl;
        Info << "Wave theory: " << waveTheory_ << endl;
        Info << "H: " << waveHeight_ << endl;
        Info << "h: " << RealwaterDepth_ << endl;
        Info << "Direction: " << waveDir_ << "º" <<  "\n\n" << endl;
    }
    else
    {
        FatalError << "Wave theory not supported, use:\n" 
	    << "Boussinesq" << exit(FatalError);
    }
