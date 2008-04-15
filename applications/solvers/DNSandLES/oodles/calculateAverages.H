    nAveragingSteps++;

    scalar nm1Coeff = scalar(nAveragingSteps - 1)/nAveragingSteps;
    scalar nCoeff = 1.0/nAveragingSteps;

    R += sqr(Umean);
    Umean = nm1Coeff*Umean + nCoeff*U;
    R = nm1Coeff*R + nCoeff*sqr(U) - sqr(Umean);
    Bmean = nm1Coeff*Bmean + nCoeff*sgsModel->B();
    epsilonMean = nm1Coeff*epsilonMean + nCoeff*sgsModel->epsilon();

    pPrime2Mean += sqr(pMean);
    pMean = nm1Coeff*pMean + nCoeff*p;
    pPrime2Mean = nm1Coeff*pPrime2Mean + nCoeff*sqr(p) - sqr(pMean);
