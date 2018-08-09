Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

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

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


//- Overset specific

// Add solver-specific interpolations
{
    dictionary suppressDict;
    {
        const wordHashSet& nonInt = Stencil::New(mesh).nonInterpolatedFields();
        for (auto fldName : nonInt)
        {
            suppressDict.add(fldName, true);
        }
        suppressDict.add("HbyA", true);
        suppressDict.add("grad(p)", true);
        suppressDict.add("surfaceIntegrate(phi)", true);
        suppressDict.add("surfaceIntegrate(phiHbyA)", true);
        suppressDict.add("cellMask", true);
        suppressDict.add("cellDisplacement", true);
        suppressDict.add("interpolatedCells", true);
        suppressDict.add("cellInterpolationWeight", true);
    }

    const_cast<dictionary&>
    (
        mesh.schemesDict()
    ).add
    (
        "oversetInterpolationSuppressed",
        suppressDict,
        true
    );
}

// Mask field for zeroing out contributions on hole cells
#include "createCellMask.H"

// Create bool field with interpolated cells
#include "createInterpolatedCells.H"

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

surfaceScalarField phiHbyA
(
    IOobject
    (
        "phiHbyA",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(U)
);

