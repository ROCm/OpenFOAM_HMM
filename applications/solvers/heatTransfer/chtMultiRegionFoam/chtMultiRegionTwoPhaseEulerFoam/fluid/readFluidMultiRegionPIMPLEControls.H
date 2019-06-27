    const dictionary& pimpleDict = mesh.solutionDict().subDict("PIMPLE");

    Switch faceMomentum
    (
        pimpleDict.lookupOrDefault<Switch>("faceMomentum", false)
    );

    int nEnergyCorrectors
    (
        pimpleDict.lookupOrDefault<int>("nEnergyCorrectors", 1)
    );
