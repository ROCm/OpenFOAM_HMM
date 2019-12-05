/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "sampledTriSurfaceMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledTriSurfaceMesh::sampleOnFaces
(
    const interpolation<Type>& sampler
) const
{
    const labelList& elements = sampleElements_;

    if (!onBoundary())
    {
        // Sample cells

        return sampledSurface::sampleOnFaces
        (
            sampler,
            elements,
            faces(),
            points()
        );
    }


    //
    // Sample boundary faces
    //

    auto tvalues = tmp<Field<Type>>::New(elements.size());
    auto& values = tvalues.ref();

    const polyBoundaryMesh& pbm = mesh().boundaryMesh();
    const label nBnd = mesh().nBoundaryFaces();

    // Create flat boundary field

    Field<Type> bVals(nBnd, Zero);

    const auto& bField = sampler.psi().boundaryField();

    forAll(bField, patchi)
    {
        const label bFacei = (pbm[patchi].start() - mesh().nInternalFaces());

        SubList<Type>
        (
            bVals,
            bField[patchi].size(),
            bFacei
        ) = bField[patchi];
    }

    // Sample in flat boundary field

    forAll(elements, i)
    {
        const label bFacei = (elements[i] - mesh().nInternalFaces());
        values[i] = bVals[bFacei];
    }

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledTriSurfaceMesh::sampleOnPoints
(
    const interpolation<Type>& interpolator
) const
{
    // One value per vertex
    auto tvalues = tmp<Field<Type>>::New(sampleElements_.size());
    auto& values = tvalues.ref();

    if (!onBoundary())
    {
        // Sample cells

        forAll(sampleElements_, pointi)
        {
            values[pointi] = interpolator.interpolate
            (
                samplePoints_[pointi],
                sampleElements_[pointi]
            );
        }
    }
    else
    {
        // Sample boundary faces

        forAll(samplePoints_, pointi)
        {
            const label facei = sampleElements_[pointi];

            values[pointi] = interpolator.interpolate
            (
                samplePoints_[pointi],
                mesh().faceOwner()[facei],
                facei
            );
        }
    }

    return tvalues;
}


// ************************************************************************* //
