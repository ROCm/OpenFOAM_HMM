Info << "Calculating averages" << endl;

scalar alpha = 
     (runTime.value() - timeToStartAveraging - runTime.deltaT().value())
    /(runTime.value() - timeToStartAveraging);

scalar onemAlpha = 1.0 - alpha;

Umean == alpha*Umean + onemAlpha*U;
rhoMean == alpha*rhoMean + onemAlpha*rho;
pMean == alpha*pMean + onemAlpha*p;
gammaMean == alpha*gammaMean + onemAlpha*gamma;
