/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "interpolatedPatch.H"
#include "meshSearch.H"
#include "polyMesh.H"
#include "interpolation.H"
#include "dictionary.H"
#include "polyPatch.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(interpolatedPatch, 0);

addToRunTimeSelectionTable(surface, interpolatedPatch, word);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class Type>
Foam::tmp<Foam::Field<Type> > Foam::interpolatedPatch::interpolate
(
    const word& fieldName,
    const fieldsCache<Type>& cache,
    const volPointInterpolation& pInterp,
    const dictionary& interpolationSchemes
) const
{
    // One value per vertex
    tmp<Field<Type> > tvalues(new Field<Type>(points().size()));
    Field<Type>& values = tvalues();

    // Get interpolator from cache
    const interpolation<Type>& interpolator =
        cache.interpolator(fieldName, pInterp, interpolationSchemes);

    if (patchIndex() != -1)
    {
        const polyPatch& patch = mesh().boundaryMesh()[patchIndex()];
        const labelList& own = mesh().faceOwner();

        boolList pointDone(points().size(), false);

        forAll(faces(), cutFaceI)
        {
            const face& f = faces()[cutFaceI];

            forAll(f, faceVertI)
            {
                label pointI = f[faceVertI];

                if (!pointDone[pointI])
                {
                    label faceI = patchFaceLabels()[cutFaceI] + patch.start();

                    label cellI = own[faceI];

                    values[pointI] =
                        interpolator.interpolate
                        (
                            points()[pointI],
                            cellI,
                            faceI
                        );
                    pointDone[pointI] = true;
                }
            }
        }
    }

    return tvalues;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interpolatedPatch::interpolatedPatch
(
    const polyMesh& mesh,
    meshSearch& searchEngine,
    const word& name,
    const word& patchName,
    const bool triangulate
)
:
    constantPatch(mesh, searchEngine, name, patchName, triangulate)
{}


Foam::interpolatedPatch::interpolatedPatch
(
    const polyMesh& mesh,
    meshSearch& searchEngine,
    const dictionary& dict          
)
:
    constantPatch(mesh, searchEngine, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interpolatedPatch::~interpolatedPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::interpolatedPatch::interpolate
(
    const word& fieldName,
    const fieldsCache<scalar>& cache,
    const volPointInterpolation& pInterp,
    const dictionary& interpolationSchemes
) const
{
    return interpolate<scalar>(fieldName, cache, pInterp, interpolationSchemes);
}


Foam::tmp<Foam::vectorField> Foam::interpolatedPatch::interpolate
(
    const word& fieldName,
    const fieldsCache<vector>& cache,
    const volPointInterpolation& pInterp,
    const dictionary& interpolationSchemes
) const
{
    return interpolate<vector>(fieldName, cache, pInterp, interpolationSchemes);
}


Foam::tmp<Foam::sphericalTensorField> Foam::interpolatedPatch::interpolate
(
    const word& fieldName,
    const fieldsCache<sphericalTensor>& cache,
    const volPointInterpolation& pInterp,
    const dictionary& interpolationSchemes
) const
{
    return interpolate<sphericalTensor>
    (
        fieldName,
        cache,
        pInterp,
        interpolationSchemes
    );
}


Foam::tmp<Foam::symmTensorField> Foam::interpolatedPatch::interpolate
(
    const word& fieldName,
    const fieldsCache<symmTensor>& cache,
    const volPointInterpolation& pInterp,
    const dictionary& interpolationSchemes
) const
{
    return interpolate<symmTensor>
    (
        fieldName,
        cache,
        pInterp,
        interpolationSchemes
    );
}


Foam::tmp<Foam::tensorField> Foam::interpolatedPatch::interpolate
(
    const word& fieldName,
    const fieldsCache<tensor>& cache,
    const volPointInterpolation& pInterp,
    const dictionary& interpolationSchemes
) const
{
    return interpolate<tensor>(fieldName, cache, pInterp, interpolationSchemes);
}


// ************************************************************************* //
