    nAveragingSteps++;

    scalar nm1Coeff = scalar(nAveragingSteps - 1)/nAveragingSteps;
    scalar nCoeff = 1.0/nAveragingSteps;

    R += sqr(Umean);
    Umean = nm1Coeff*Umean + nCoeff*U;
    R = nm1Coeff*R + nCoeff*sqr(U) - sqr(Umean);

    pPrime2Mean += sqr(pMean);
    pMean = nm1Coeff*pMean + nCoeff*p;
    pPrime2Mean = nm1Coeff*pPrime2Mean + nCoeff*sqr(p) - sqr(pMean);

    Tmean = nm1Coeff*Tmean + nCoeff*T;
    bmean = nm1Coeff*bmean + nCoeff*b;
    Ximean = nm1Coeff*Ximean + nCoeff*Xi;
