    const dictionary& simple = mesh.solutionDict().subDict("SIMPLE");

    const int nNonOrthCorr =
        simple.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

    const bool momentumPredictor =
        simple.lookupOrDefault("momentumPredictor", true);

    const bool transonic =
        simple.lookupOrDefault("transonic", false);

