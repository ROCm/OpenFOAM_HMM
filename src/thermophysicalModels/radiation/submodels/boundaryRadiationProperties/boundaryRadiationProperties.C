/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenCFD Ltd
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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::radiation::boundaryRadiationProperties::createIOobject
(
    const fvMesh& mesh, const word name
) const
{
    IOobject io
    (
        name,
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.headerOk())
    {
        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
        return io;
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
    const IOobject boundaryIO
    (
        createIOobject(mesh, boundaryRadiationProperties::typeName)
    );

    if (boundaryIO.readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        radBoundaryProperties_.set
        (
            new volScalarField(boundaryIO, mesh)
        );
    }
}


// * * * * * * * * * * * * * * * Member fucntions * * * * * * * * * * * * *  //

const Foam::volScalarField& Foam::radiation::boundaryRadiationProperties::
radBoundaryProperties() const
{
    return radBoundaryProperties_();
}


Foam::tmp<Foam::scalarField> Foam::radiation::boundaryRadiationProperties::
emissivity(const label index, const label bandI) const
{
    if (!radBoundaryProperties_.empty())
    {
        return refCast<const boundaryRadiationPropertiesFvPatchField>
        (
            radBoundaryProperties_->boundaryField()[index]
        ).emissivity(bandI);
    }
    else
    {
        FatalErrorInFunction
            << "Field 'boundaryRadiationProperties'"
            << "is not found in the constant directory."
            << "Please add it "
            << exit(FatalError);

         return tmp<scalarField>(new scalarField());
    }
}


Foam::tmp<Foam::scalarField> Foam::radiation::boundaryRadiationProperties::
absorptivity(const label index, const label bandI) const
{
    if (!radBoundaryProperties_.empty())
    {
        return refCast<const boundaryRadiationPropertiesFvPatchField>
        (
            radBoundaryProperties_->boundaryField()[index]
        ).absorptivity(bandI);
    }
    else
    {
        FatalErrorInFunction
            << "Field 'boundaryRadiationProperties'"
            << "is not found in the constant directory."
            << "Please add it "
            << exit(FatalError);

        return tmp<scalarField>(new scalarField());
    }
}


Foam::tmp<Foam::scalarField> Foam::radiation::boundaryRadiationProperties::
transmissivity(const label index, const label bandI) const
{
    if (!radBoundaryProperties_.empty())
    {
        return refCast<const boundaryRadiationPropertiesFvPatchField>
        (
            radBoundaryProperties_->boundaryField()[index]
        ).transmissivity(bandI);
    }
    else
    {
        FatalErrorInFunction
            << "Field 'boundaryRadiationProperties'"
            << "is not found in the constant directory."
            << "Please add it "
            << exit(FatalError);

        return tmp<scalarField>(new scalarField());
    }
}


Foam::tmp<Foam::scalarField> Foam::radiation::boundaryRadiationProperties::
reflectivity(const label index, const label bandI) const
{
    if (!radBoundaryProperties_.empty())
    {
        return refCast<const boundaryRadiationPropertiesFvPatchField>
        (
            radBoundaryProperties_->boundaryField()[index]
        ).reflectivity(bandI);
    }
    else
    {
        FatalErrorInFunction
            << "Field 'boundaryRadiationProperties'"
            << "is not found in the constant directory."
            << "Please add it "
            << exit(FatalError);

        return tmp<scalarField>(new scalarField());
    }
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::radiation::boundaryRadiationProperties::~boundaryRadiationProperties()
{}


// ************************************************************************* //
