    nAveragingSteps++;

    scalar nm1Coeff = scalar(nAveragingSteps - 1)/nAveragingSteps;
    scalar nCoeff = 1.0/nAveragingSteps;

    //R += sqr(Umean);
    kGSMean += 0.5*tr(sqr(Umean));
    Umean = nm1Coeff*Umean + nCoeff*U;
    kGSMean = nm1Coeff*kGSMean + nCoeff*0.5*tr(sqr(U)) - 0.5*tr(sqr(Umean));
    //R = nm1Coeff*R + nCoeff*sqr(U) - sqr(Umean);
    //Bmean = nm1Coeff*Bmean + nCoeff*turbulence->B();
    kMean = nm1Coeff*kMean + nCoeff*turbulence->k();
    nuEffMean = nm1Coeff*nuEffMean + nCoeff*turbulence->nuEff();
    gammaMean = nm1Coeff*gammaMean + nCoeff*gamma;
    sigmaMean = nm1Coeff*sigmaMean + nCoeff*mag(fvc::grad(gamma));

    pPrime2Mean += sqr(pMean);
    pMean = nm1Coeff*pMean + nCoeff*pd;
    pPrime2Mean = nm1Coeff*pPrime2Mean + nCoeff*sqr(pd) - sqr(pMean);
