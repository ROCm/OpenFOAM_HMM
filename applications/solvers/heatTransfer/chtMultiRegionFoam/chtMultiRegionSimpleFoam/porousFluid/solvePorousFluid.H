//  Pressure-velocity SIMPLE corrector

    porousP.storePrevIter();
    porousRho.storePrevIter();
    {
        #include "UPorousFluidEqn.H"
        #include "hPorousFluidEqn.H"
        #include "pPorousFluidEqn.H"
    }

    turbPorous.correct();
