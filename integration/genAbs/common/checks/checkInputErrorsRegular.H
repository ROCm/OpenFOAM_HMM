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
else if ( waveTheory_ == "StokesI" || waveTheory_ == "StokesII" )
{
    if ( wavePeriod_ <= 0.0 )
    {
        FatalError << "Check wave period value." << exit(FatalError);
    }

    if ( waveLength_ <= 0.0 )
    {
        word generation = "a";

        waveLength_ = StokesIFun::waveLength (RealwaterDepth_, wavePeriod_);

        waveK = 2.0*PI()/waveLength_;

        if ( waveK*RealwaterDepth_ > PI() )
        {
            generation = "Deep";
        }
        else if ( waveK*RealwaterDepth_ < PI()/10.0 )
        {
            generation = "Shallow";
        }
        else
        {
            generation = "Intermediate";
        }

        Info << "\nIH Wave Generation BC" << endl;
        Info << "	Wave theory: " << waveTheory_ << endl;
        Info << "	H: " << waveHeight_ << endl;
        Info << "	T: " << wavePeriod_ << endl;
        Info << "	h: " << RealwaterDepth_ << endl;
        Info << "	waveLength: " << waveLength_ << endl;
        Info << "	Direction: " << waveDir_ << "º" << endl;
        Info << "	Generation in: " << generation << " waters." << endl;
        Info << "	Relative depth (kh): " << waveK*RealwaterDepth_ 
	     << "\n\n" << endl;
    }
}
else if ( waveTheory_ == "StokesV" )
{
    if ( wavePeriod_ <= 0.0 )
    {
        FatalError << "Check wave period value." << exit(FatalError);
    }

    if ( lambdaStokesV_ <= 0.0 )
    {
        word generation = "a";

        scalar f1;
        scalar f2;

        stokesVFun :: StokesVNR ( waveHeight_, 
            RealwaterDepth_,  
            wavePeriod_, 
            &waveK, 
            &lambdaStokesV_, 
            &f1, 
            &f2 
            );

        Info << "f1 residual " << f1 << endl;
        Info << "f2 residual " << f2 << endl;

        if ( f1 > 0.001 || f2 > 0.001 )
        {
	    FatalError << "No convergence for Stokes V wave theory.\n" 
	        << exit(FatalError);
        }

        waveLength_ = 2.0*PI()/waveK;

        if ( waveK*RealwaterDepth_ > PI() )
        {
            generation = "Deep";
        }
        else if ( waveK*RealwaterDepth_ < PI()/10.0 )
        {
            generation = "Shallow";
        }
        else
        {
            generation = "Intermediate";
        }

        Info << "\nIH Wave Generation BC" << endl;
        Info << "	Wave theory: " << waveTheory_ << endl;
        Info << "	H: " << waveHeight_ << endl;
        Info << "	T: " << wavePeriod_ << endl;
        Info << "	h: " << RealwaterDepth_ << endl;
        Info << "	L: " << waveLength_ << endl;
        Info << "	Lambda: " << lambdaStokesV_ << endl;
        Info << "	Direction: " << waveDir_ << "º"  << endl;
        Info << "	Generation in: " << generation << " waters." << endl;
        Info << "	Relative depth (kh): " << waveK*RealwaterDepth_ 
	     << "\n\n" << endl;
    }
}
else if ( waveTheory_ == "Cnoidal" )
{
    if ( wavePeriod_ <= 0.0 )
    {
        FatalError << "Check wave period value." << exit(FatalError);
    }

    if ( mCnoidal_ <= 0.0 )
    {
        cnoidalFun :: calculations ( waveHeight_, 
            RealwaterDepth_,  
            wavePeriod_, 
            &mCnoidal_, 
            &waveLength_);

        Info << "\nIH Wave Generation BC" << endl;
        Info << "	Wave theory: " << waveTheory_ << endl;
        Info << "	H: " << waveHeight_ << endl;
        Info << "	T: " << wavePeriod_ << endl;
        Info << "	h: " << RealwaterDepth_ << endl;
        Info << "	L: " << waveLength_ << endl;
        Info << "	m parameter: " << mCnoidal_ << endl;
        Info << "	Direction: " << waveDir_ << "º" <<  "\n\n" << endl;
    }
}
else
{
    FatalError << "Wave theory not supported, use:\n" 
        << "StokesI, StokesII, StokesV, Cnoidal." << exit(FatalError);
}

