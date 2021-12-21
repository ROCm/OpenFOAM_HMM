/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "wallHeatFlux.H"
#include "turbulentFluidThermoModel.H"
#include "solidThermo.H"
#include "surfaceInterpolate.H"
#include "fvcSnGrad.H"
#include "wallPolyPatch.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(wallHeatFlux, 0);
    addToRunTimeSelectionTable(functionObject, wallHeatFlux, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::wallHeatFlux::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Wall heat-flux");
    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    writeTabbed(os, "integral");
    os  << endl;
}


void Foam::functionObjects::wallHeatFlux::calcHeatFlux
(
    const volScalarField& alpha,
    const volScalarField& he,
    volScalarField& wallHeatFlux
)
{
    volScalarField::Boundary& wallHeatFluxBf = wallHeatFlux.boundaryFieldRef();

    const volScalarField::Boundary& heBf = he.boundaryField();

    const volScalarField::Boundary& alphaBf = alpha.boundaryField();

    for (const label patchi : patchSet_)
    {
        wallHeatFluxBf[patchi] = alphaBf[patchi]*heBf[patchi].snGrad();
    }


    const auto* qrPtr = cfindObject<volScalarField>(qrName_);

    if (qrPtr)
    {
        const volScalarField::Boundary& radHeatFluxBf = qrPtr->boundaryField();

        for (const label patchi : patchSet_)
        {
            wallHeatFluxBf[patchi] -= radHeatFluxBf[patchi];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::wallHeatFlux::wallHeatFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    patchSet_(),
    qrName_("qr")
{
    volScalarField* wallHeatFluxPtr
    (
        new volScalarField
        (
            IOobject
            (
                scopedName(typeName),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimMass/pow3(dimTime), Zero)
        )
    );

    mesh_.objectRegistry::store(wallHeatFluxPtr);

    read(dict);

    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::wallHeatFlux::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            dict.getOrDefault<wordRes>("patches", wordRes())
        );

    dict.readIfPresent("qr", qrName_);

    Info<< type() << " " << name() << ":" << nl;

    if (patchSet_.empty())
    {
        forAll(pbm, patchi)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                patchSet_.insert(patchi);
            }
        }

        Info<< "    processing all wall patches" << nl << endl;
    }
    else
    {
        Info<< "    processing wall patches: " << nl;
        labelHashSet filteredPatchSet;
        for (const label patchi : patchSet_)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                filteredPatchSet.insert(patchi);
                Info<< "        " << pbm[patchi].name() << endl;
            }
            else
            {
                WarningInFunction
                    << "Requested wall heat-flux on non-wall boundary "
                    << "type patch: " << pbm[patchi].name() << endl;
            }
        }

        Info<< endl;

        patchSet_ = filteredPatchSet;
    }

    return true;
}


bool Foam::functionObjects::wallHeatFlux::execute()
{
    auto& wallHeatFlux = lookupObjectRef<volScalarField>(scopedName(typeName));

    if
    (
        foundObject<compressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        )
    )
    {
        const compressible::turbulenceModel& turbModel =
            lookupObject<compressible::turbulenceModel>
            (
                turbulenceModel::propertiesName
            );

        calcHeatFlux
        (
            turbModel.alphaEff()(),
            turbModel.transport().he(),
            wallHeatFlux
        );
    }
    else if (foundObject<fluidThermo>(fluidThermo::dictName))
    {
        const fluidThermo& thermo =
            lookupObject<fluidThermo>(fluidThermo::dictName);

        calcHeatFlux
        (
            thermo.alpha(),
            thermo.he(),
            wallHeatFlux
        );
    }
    else if (foundObject<solidThermo>(solidThermo::dictName))
    {
        const solidThermo& thermo =
            lookupObject<solidThermo>(solidThermo::dictName);

        calcHeatFlux(thermo.alpha(), thermo.he(), wallHeatFlux);
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find compressible turbulence model in the "
            << "database" << exit(FatalError);
    }

    const fvPatchList& patches = mesh_.boundary();

    const surfaceScalarField::Boundary& magSf = mesh_.magSf().boundaryField();

    for (const label patchi : patchSet_)
    {
        const fvPatch& pp = patches[patchi];

        const scalarField& hfp = wallHeatFlux.boundaryField()[patchi];

        const scalar minHfp = gMin(hfp);
        const scalar maxHfp = gMax(hfp);
        const scalar integralHfp = gSum(magSf[patchi]*hfp);

        if (Pstream::master())
        {
            writeCurrentTime(file());

            file()
                << token::TAB << pp.name()
                << token::TAB << minHfp
                << token::TAB << maxHfp
                << token::TAB << integralHfp
                << endl;
        }

        Log << "    min/max/integ(" << pp.name() << ") = "
            << minHfp << ", " << maxHfp << ", " << integralHfp << endl;

        this->setResult("min(" + pp.name() + ")", minHfp);
        this->setResult("max(" + pp.name() + ")", maxHfp);
        this->setResult("int(" + pp.name() + ")", integralHfp);
    }


    return true;
}


bool Foam::functionObjects::wallHeatFlux::write()
{
    const auto& wallHeatFlux =
        lookupObject<volScalarField>(scopedName(typeName));

    Log << type() << " " << name() << " write:" << nl
        << "    writing field " << wallHeatFlux.name() << endl;

    wallHeatFlux.write();

    return true;
}


// ************************************************************************* //
