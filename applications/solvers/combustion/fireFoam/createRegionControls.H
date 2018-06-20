IOdictionary additionalControlsDict
(
    IOobject
    (
        "additionalControls",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

bool solvePrimaryRegion
(
    additionalControlsDict.lookupOrDefault("solvePrimaryRegion", true)
);

bool solvePyrolysisRegion
(
    additionalControlsDict.lookupOrDefault("solvePyrolysisRegion", true)
);

scalar maxDi = pyrolysis.maxDiff();
