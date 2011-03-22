    // We do not have a top-level mesh. Construct the fvSolution for
    // the runTime instead.
    runTime.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
    fvSolution solutionDict(runTime);

    const dictionary& pimple = solutionDict.subDict("PIMPLE");

    const int nOuterCorr =
        pimple.lookupOrDefault<int>("nOuterCorrectors", 1);
