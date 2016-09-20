/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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

#include "PecletNo.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "surfaceInterpolate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(PecletNo, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        PecletNo,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::PecletNo::calc()
{
    if (foundObject<surfaceScalarField>(fieldName_))
    {
        // Obtain nuEff of muEff.  Assumes that a compressible flux is present
        // when using a compressible turbulence model, and an incompressible
        // flux when using an incompressible turbulence model

        tmp<volScalarField> nuOrMuEff;
        if (mesh_.foundObject<cmpTurbModel>(turbulenceModel::propertiesName))
        {
            const cmpTurbModel& model =
                mesh_.lookupObject<cmpTurbModel>
                (
                    turbulenceModel::propertiesName
                );

            nuOrMuEff = model.muEff();
        }
        else if
        (
            mesh_.foundObject<icoTurbModel>(turbulenceModel::propertiesName)
        )
        {
            const icoTurbModel& model =
                mesh_.lookupObject<icoTurbModel>
                (
                    turbulenceModel::propertiesName
                );

            nuOrMuEff = model.nuEff();
        }
        else if (mesh_.foundObject<dictionary>("transportProperties"))
        {
            const dictionary& model =
                mesh_.lookupObject<dictionary>("transportProperties");

            nuOrMuEff =
                tmp<volScalarField>
                (
                    new volScalarField
                    (
                        IOobject
                        (
                            "nuEff",
                            mesh_.time().timeName(),
                            mesh_,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh_,
                        dimensionedScalar(model.lookup("nu"))
                    )
                );
        }
        else
        {
            FatalErrorInFunction
                << "Unable to determine the viscosity"
                << exit(FatalError);
        }


        const surfaceScalarField& phi =
            mesh_.lookupObject<surfaceScalarField>(fieldName_);

        return store
        (
            resultName_,
            mag(phi)
           /(
                mesh_.magSf()
               *mesh_.surfaceInterpolation::deltaCoeffs()
               *fvc::interpolate(nuOrMuEff)
            )
        );
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::PecletNo::PecletNo
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "phi")
{
    setResultName("Pe", "phi");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::PecletNo::~PecletNo()
{}


// ************************************************************************* //
