/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

#include "boundaryRadiationProperties.H"
#include "boundaryRadiationPropertiesFvPatchField.H"
#include "fvPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(boundaryRadiationProperties, 0);
    }
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::radiation::boundaryRadiationProperties::boundaryRadiationProperties
(
    const fvMesh& mesh
)
:
    MeshObject
    <
        fvMesh,
        Foam::GeometricMeshObject,
        boundaryRadiationProperties
    >(mesh),
    radBoundaryProperties_()
{
    IOobject boundaryIO
    (
        boundaryRadiationProperties::typeName,
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (boundaryIO.typeHeaderOk<volScalarField>(true))
    {
        radBoundaryProperties_.set
        (
            new volScalarField(boundaryIO, mesh)
        );
    }
}


// * * * * * * * * * * * * * * * Member fucntions * * * * * * * * * * * * *  //

Foam::tmp<Foam::scalarField>
Foam::radiation::boundaryRadiationProperties::emissivity
(
    const label patchI,
    const label bandI
) const
{
    if (!radBoundaryProperties_.empty())
    {
        return refCast<const boundaryRadiationPropertiesFvPatchField>
        (
            radBoundaryProperties_->boundaryField()[patchI]
        ).emissivity(bandI);
    }
    else
    {
        FatalErrorInFunction
            << "Field 'boundaryRadiationProperties'"
            << "is not found in the constant directory. "
            << "Please add it"
            << exit(FatalError);

        return tmp<scalarField>(new scalarField());
    }
}


Foam::tmp<Foam::scalarField>
Foam::radiation::boundaryRadiationProperties::absorptivity
(
    const label patchI,
    const label bandI
) const
{
    if (!radBoundaryProperties_.empty())
    {
        return refCast<const boundaryRadiationPropertiesFvPatchField>
        (
            radBoundaryProperties_->boundaryField()[patchI]
        ).absorptivity(bandI);
    }
    else
    {
        FatalErrorInFunction
            << "Field 'boundaryRadiationProperties'"
            << "is not found in the constant directory. "
            << "Please add it "
            << exit(FatalError);

        return tmp<scalarField>(new scalarField());
    }
}


Foam::tmp<Foam::scalarField>
Foam::radiation::boundaryRadiationProperties::transmissivity
(
    const label patchI,
    const label bandI
) const
{
    if (!radBoundaryProperties_.empty())
    {
        return refCast<const boundaryRadiationPropertiesFvPatchField>
        (
            radBoundaryProperties_->boundaryField()[patchI]
        ).transmissivity(bandI);
    }
    else
    {
        FatalErrorInFunction
            << "Field 'boundaryRadiationProperties'"
            << "is not found in the constant directory. "
            << "Please add it"
            << exit(FatalError);

        return tmp<scalarField>(new scalarField());
    }
}


Foam::tmp<Foam::scalarField>
Foam::radiation::boundaryRadiationProperties::reflectivity
(
    const label patchI,
    const label bandI
) const
{
    if (!radBoundaryProperties_.empty())
    {
        return refCast<const boundaryRadiationPropertiesFvPatchField>
        (
            radBoundaryProperties_->boundaryField()[patchI]
        ).reflectivity(bandI);
    }
    else
    {
        FatalErrorInFunction
            << "Field 'boundaryRadiationProperties'"
            << "is not found in the constant directory. "
            << "Please add it"
            << exit(FatalError);

        return tmp<scalarField>(new scalarField());
    }
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::radiation::boundaryRadiationProperties::~boundaryRadiationProperties()
{}


// ************************************************************************* //
