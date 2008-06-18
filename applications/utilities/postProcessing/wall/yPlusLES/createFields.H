Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

#include "createPhi.H"

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::LESmodel> sgsModel
(
    incompressible::LESmodel::New(U, phi, laminarTransport)
);

volScalarField::GeometricBoundaryField d = nearWallDist(mesh).y();
