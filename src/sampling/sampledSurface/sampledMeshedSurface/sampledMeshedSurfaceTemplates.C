/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "sampledMeshedSurface.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledMeshedSurface::sampleOnFaces
(
    const interpolation<Type>& sampler
) const
{
    const Type deflt
    (
        defaultValues_.getOrDefault<Type>(sampler.psi().name(), Zero)
    );

    const labelList& elements = sampleElements_;


    if (!onBoundary())
    {
        // Sample cells

        return sampledSurface::sampleOnFaces
        (
            sampler,
            elements,
            faces(),
            points(),
            deflt
        );
    }


    //
    // Sample boundary faces
    //

    auto tvalues = tmp<Field<Type>>::New(elements.size());
    auto& values = tvalues.ref();

    // Create flat boundary field

    const polyBoundaryMesh& pbm = mesh().boundaryMesh();

    Field<Type> bVals(mesh().nBoundaryFaces(), Zero);

    const auto& bField = sampler.psi().boundaryField();

    forAll(bField, patchi)
    {
        SubList<Type>(bVals, pbm[patchi].range()) = bField[patchi];
    }

    // Sample within the flat boundary field

    forAll(elements, i)
    {
        const label bFacei = (elements[i] - mesh().nInternalFaces());

        if (bFacei < 0)
        {
            values[i] = deflt;
        }
        else
        {
            values[i] = bVals[bFacei];
        }
    }

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledMeshedSurface::sampleOnPoints
(
    const interpolation<Type>& interpolator
) const
{
    const Type deflt
    (
        defaultValues_.getOrDefault<Type>(interpolator.psi().name(), Zero)
    );

    const labelList& elements = sampleElements_;


    // One value per vertex.
    // - sampleElements_ and samplePoints_ have the identical size
    auto tvalues = tmp<Field<Type>>::New(elements.size());
    auto& values = tvalues.ref();

    if (!onBoundary())
    {
        // Sample cells

        forAll(elements, i)
        {
            const label celli = elements[i];

            if (celli < 0)
            {
                values[i] = deflt;
            }
            else
            {
                values[i] = interpolator.interpolate
                (
                    samplePoints_[i],
                    celli
                );
            }
        }

        return tvalues;
    }


    //
    // Sample boundary faces
    //

    forAll(elements, i)
    {
        const label facei = elements[i];

        if (facei < 0)
        {
            values[i] = deflt;
        }
        else
        {
            values[i] = interpolator.interpolate
            (
                samplePoints_[i],
                mesh().faceOwner()[facei],
                facei
            );
        }
    }

    return tvalues;
}


// ************************************************************************* //
