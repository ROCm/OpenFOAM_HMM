/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "MovingPhaseModel.H"

#include "phaseSystem.H"

#include "fixedValueFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"

#include "fvmDdt.H"
#include "fvmDiv.H"

#include "fvmSup.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::generatePhi
(
    const word& phiName,
    const volVectorField& U
) const
{
    IOobject phiHeader
    (
        phiName,
        U.mesh().time().timeName(),
        U.mesh(),
        IOobject::NO_READ
    );

    if (phiHeader.typeHeaderOk<surfaceScalarField>(true))
    {
        Info<< "Reading face flux field " << phiName << endl;

        return
            tmp<surfaceScalarField>
            (
                new surfaceScalarField
                (
                    IOobject
                    (
                        phiName,
                        U.mesh().time().timeName(),
                        U.mesh(),
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    U.mesh()
                )
            );
    }
    else
    {
        Info<< "Calculating face flux field " << phiName << endl;

        wordList phiTypes
        (
            U.boundaryField().size(),
            calculatedFvPatchScalarField::typeName
        );

        forAll(U.boundaryField(), i)
        {
            if
            (
                isA<fixedValueFvPatchVectorField>(U.boundaryField()[i])
             || isA<slipFvPatchVectorField>(U.boundaryField()[i])
             || isA<partialSlipFvPatchVectorField>(U.boundaryField()[i])
            )
            {
                phiTypes[i] = fixedValueFvPatchScalarField::typeName;
            }
        }

        return tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    U.mesh().time().timeName(),
                    U.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fvc::interpolate(U) & U.mesh().Sf(),
                phiTypes
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::MovingPhaseModel<BasePhaseModel>::MovingPhaseModel
(
    const phaseSystem& fluid,
    const word& phaseName
)
:
    BasePhaseModel(fluid, phaseName),
    U_(fluid.mesh().lookupObject<volVectorField>("U")),
    phi_(fluid.mesh().lookupObject<surfaceScalarField>("phi")),
    /*
    (
        generatePhi
        (
            IOobject::groupName("phi", phaseModel::name()),
            U_
        )
    ),
    */
    alphaPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaPhi", phaseModel::name()),
            fluid.mesh().time().timeName(),
            fluid.mesh()
        ),
        fluid.mesh(),
        dimensionedScalar("0", dimensionSet(0, 3, -1, 0, 0), 0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseModel>
Foam::MovingPhaseModel<BasePhaseModel>::~MovingPhaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasePhaseModel>
void Foam::MovingPhaseModel<BasePhaseModel>::correct()
{
    BasePhaseModel::correct();
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::phi() const
{
    return tmp<surfaceScalarField>(phi_);
}


template<class BasePhaseModel>
const Foam::surfaceScalarField&
Foam::MovingPhaseModel<BasePhaseModel>::phi()
{
    return phi_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::MovingPhaseModel<BasePhaseModel>::alphaPhi() const
{
    return tmp<surfaceScalarField>(alphaPhi_);
}


template<class BasePhaseModel>
Foam::surfaceScalarField&
Foam::MovingPhaseModel<BasePhaseModel>::alphaPhi()
{
    return alphaPhi_;
}


template<class BasePhaseModel>
Foam::tmp<Foam::volVectorField>
Foam::MovingPhaseModel<BasePhaseModel>::U() const
{
    return tmp<volVectorField>(U_);
}


// ************************************************************************* //
