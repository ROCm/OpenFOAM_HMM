{

    // calculate the new diameter with the standard 1D Rosin Rammler distribution
    //-------------AL_____101012---------------------------//
    // Calculation of the mean radius based on SMR rs. Coefficient factorGamma depends on nExp.
    // Note that Reitz either used (Schmidt et al., 1999-01-0496) or skipped (Senecal et al.) this factor!!!
    //      scalar factorGamma = 0.75*sqrt(mathematicalConstant::pi);       //nExp=2
    scalar factorGamma = 1.;
    scalar delta = dD/factorGamma;

    /*      dD is the SMD, and the delta is calculated using gama
        function. Here we assume nExp = 2. */
    //      scalar delta = dD/(0.75*sqrt(mathematicalConstant::pi));

    //        scalar minValue = min(p.d()/20.0,dD/20.0);
    scalar minValue = dD/10.0;

    //      delta is divided by 20 instead of 10 in order to make sure of small minValue
    //      scalar minValue = min(p.d(),dD/20.0);

    //        scalar maxValue = p.d();
    scalar maxValue = dD;

    //      The pdf value for 4.0*delta is already very small.
    //      scalar maxValue = delta*4.0;

    if(maxValue - minValue < SMALL)
    {
        //          minValue = p.d()/20.0;
        minValue = maxValue/20.0;
        //-----------------------------------END-------------------------------------//
    }

    scalar range = maxValue - minValue;

    scalar nExp = 3;
    scalar rrd[500];
    //--------------------------------AL_____101012------------------------------//
    scalar probFactorMin = exp(-pow(minValue/delta,nExp));
    scalar probFactorMax = exp(-pow(maxValue/delta,nExp));
    scalar probFactor = 1./(probFactorMin - probFactorMax);
    //-----------------------------------END-------------------------------------//
    for(label n=0; n<500; n++)
    {
        scalar xx = minValue + range*n/500;
        //-------------------------------AL_____101012-------------------------------//
        //          rrd_[n] = 1 - exp(-pow(xx/delta,nExp));
        rrd[n] = (probFactorMin - exp(-pow(xx/delta,nExp)))*probFactor;
        //-----------------------------------END-------------------------------------//
    }


    bool success = false;

    scalar y = rndGen.scalar01();
    label k = 0;

    while(!success && (k<500))
    {
        if (rrd[k]>y)
        {
            success = true;
        }
        k++;

    }
    //--------------------------------AL_____101012------------------------------//
    //      x = minValue + range*n/500;
    x = minValue + range*(k-0.5)/500.0;
    //------------------------------------END------------------------------------//
}
