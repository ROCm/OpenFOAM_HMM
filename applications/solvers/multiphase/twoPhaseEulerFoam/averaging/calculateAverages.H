    nAveragingSteps++;

    scalar nm1Coeff = scalar(nAveragingSteps - 1)/nAveragingSteps;
    scalar nCoeff = 1.0/nAveragingSteps;

    UaMean = nm1Coeff*UaMean + nCoeff*Ua;
    UbMean = nm1Coeff*UbMean + nCoeff*Ub;

    alphaMean = nm1Coeff*alphaMean + nCoeff*alpha;
    pMean = nm1Coeff*pMean + nCoeff*p;
