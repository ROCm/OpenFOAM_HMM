/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenCFD Ltd.
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
    radBoundaryPropertiesPtrList_(mesh.boundary().size())
{
    IOobject boundaryIO
    (
        boundaryRadiationProperties::typeName,
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (boundaryIO.typeHeaderOk<IOdictionary>(true))
    {
        const IOdictionary radiationDict(boundaryIO);

        forAll(mesh.boundary(), patchi)
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchi];

            if (radiationDict.isDict(pp.name()))
            {
                const dictionary& dict = radiationDict.subDict(pp.name());

                radBoundaryPropertiesPtrList_[patchi].reset
                (
                    new boundaryRadiationPropertiesPatch(pp, dict)
                );
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * *  //

Foam::tmp<Foam::scalarField>
Foam::radiation::boundaryRadiationProperties::emissivity
(
    const label patchi,
    const label bandi
) const
{
    if (!radBoundaryPropertiesPtrList_[patchi].empty())
    {
        return radBoundaryPropertiesPtrList_[patchi]->emissivity(bandi);
    }

    FatalErrorInFunction
        << "Patch : " << mesh().boundaryMesh()[patchi].name()
        << " is not found in the boundaryRadiationProperties. "
        << "Please add it"
        << exit(FatalError);

    return tmp<scalarField>::New();
}


Foam::tmp<Foam::scalarField>
Foam::radiation::boundaryRadiationProperties::absorptivity
(
    const label patchi,
    const label bandi
) const
{
    if (!radBoundaryPropertiesPtrList_[patchi].empty())
    {
        return radBoundaryPropertiesPtrList_[patchi]->absorptivity(bandi);
    }

     FatalErrorInFunction
         << "Patch : " << mesh().boundaryMesh()[patchi].name()
         << " is not found in the boundaryRadiationProperties. "
         << "Please add it"
         << exit(FatalError);

    return tmp<scalarField>::New();
}


Foam::tmp<Foam::scalarField>
Foam::radiation::boundaryRadiationProperties::transmissivity
(
    const label patchi,
    const label bandi
) const
{
    if (!radBoundaryPropertiesPtrList_[patchi].empty())
    {
        return radBoundaryPropertiesPtrList_[patchi]->transmissivity(bandi);
    }

    FatalErrorInFunction
        << "Patch : " << mesh().boundaryMesh()[patchi].name()
        << " is not found in the boundaryRadiationProperties. "
        << "Please add it"
        << exit(FatalError);

    return tmp<scalarField>::New();
}


Foam::tmp<Foam::scalarField>
Foam::radiation::boundaryRadiationProperties::reflectivity
(
    const label patchi,
    const label bandi
) const
{
    if (!radBoundaryPropertiesPtrList_[patchi].empty())
    {
        return radBoundaryPropertiesPtrList_[patchi]->reflectivity(bandi);
    }

    FatalErrorInFunction
        << "Patch : " << mesh().boundaryMesh()[patchi].name()
        << " is not found in the boundaryRadiationProperties. "
        << "Please add it"
        << exit(FatalError);

    return tmp<scalarField>::New();
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::radiation::boundaryRadiationProperties::~boundaryRadiationProperties()
{}


// ************************************************************************* //
