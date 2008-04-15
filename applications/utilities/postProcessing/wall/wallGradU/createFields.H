    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

#   include "createPhi.H"

    autoPtr<transportModel> laminarTransport
    (
        transportModel::New(U, phi)
    );

    autoPtr<turbulenceModel> turbulence
    (
        turbulenceModel::New(U, phi, laminarTransport())
    );
