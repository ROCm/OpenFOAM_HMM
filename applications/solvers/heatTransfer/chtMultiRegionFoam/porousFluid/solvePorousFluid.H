if (finalIter)
{
    porousMesh.data::add("finalIteration", true);
}

if (oCorr == 0)
{
    solve(fvm::ddt(porousRho) + fvc::div(porousPhi));
}

#include "UPorousFluidEqn.H"

#include "hPorousFluidEqn.H"

// --- PISO loop
for (int corr=0; corr<nCorr; corr++)
{
    #include "pPorousFluidEqn.H"
}

turbPorous.correct();

porousRho = porousThermo.rho();

if (finalIter)
{
    porousMesh.data::remove("finalIteration");
}
