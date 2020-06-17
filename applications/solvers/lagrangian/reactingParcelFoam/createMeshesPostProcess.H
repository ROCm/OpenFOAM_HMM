#include "createMesh.H"

IOdictionary filmDict
(
    IOobject
    (
        "surfaceFilmProperties",
        runTime.constant(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    )
);

const word filmRegionName = filmDict.get<word>("region");

fvMesh filmMesh
(
    IOobject
    (
        filmRegionName,
        runTime.timeName(),
        runTime,
        IOobject::MUST_READ
    )
);
