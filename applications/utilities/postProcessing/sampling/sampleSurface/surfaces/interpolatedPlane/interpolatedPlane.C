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

#include "interpolatedPlane.H"
#include "meshSearch.H"
#include "polyMesh.H"
#include "meshSearch.H"
#include "interpolation.H"
#include "dictionary.H"
#include "plane.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(interpolatedPlane, 0);

addToRunTimeSelectionTable(surface, interpolatedPlane, word);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class Type>
Foam::tmp<Foam::Field<Type> > Foam::interpolatedPlane::interpolate
(
    const word& fieldName,
    const fieldsCache<Type>& cache,
    const volPointInterpolation& pInterp,
    const dictionary& interpolationSchemes
) const
{
    // One value per point
    tmp<Field<Type> > tvalues(new Field<Type>(points().size()));
    Field<Type>& values = tvalues();

    // Get interpolator from cache
    const interpolation<Type>& interpolator =
        cache.interpolator(fieldName, pInterp, interpolationSchemes);

    boolList pointDone(points().size(), false);

    forAll(faces(), cutFaceI)
    {
        const face& f = faces()[cutFaceI];

        forAll(f, faceVertI)
        {
            label pointI = f[faceVertI];

            if (!pointDone[pointI])
            {
                values[pointI] =
                    interpolator.interpolate
                    (
                        points()[pointI],
                        meshCells()[cutFaceI]
                    );
                pointDone[pointI] = true;
            }
        }
    }

    return tvalues;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::interpolatedPlane::interpolatedPlane
(
    const polyMesh& mesh,
    meshSearch& searchEngine,
    const word& name,
    const plane& planeDesc,
    const bool triangulate
)
:
    constantPlane(mesh, searchEngine, name, planeDesc, triangulate)
{}


// Construct from dictionary
Foam::interpolatedPlane::interpolatedPlane
(
    const polyMesh& mesh,
    meshSearch& searchEngine,
    const dictionary& dict          
)
:
    constantPlane(mesh, searchEngine, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interpolatedPlane::~interpolatedPlane()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::interpolatedPlane::interpolate
(
    const word& fieldName,
    const fieldsCache<scalar>& cache,
    const volPointInterpolation& pInterp,
    const dictionary& interpolationSchemes
) const
{
    return interpolate<scalar>(fieldName, cache, pInterp, interpolationSchemes);
}


Foam::tmp<Foam::vectorField> Foam::interpolatedPlane::interpolate
(
    const word& fieldName,
    const fieldsCache<vector>& cache,
    const volPointInterpolation& pInterp,
    const dictionary& interpolationSchemes
) const
{
    return interpolate<vector>(fieldName, cache, pInterp, interpolationSchemes);
}


Foam::tmp<Foam::sphericalTensorField> Foam::interpolatedPlane::interpolate
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


Foam::tmp<Foam::symmTensorField> Foam::interpolatedPlane::interpolate
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


Foam::tmp<Foam::tensorField> Foam::interpolatedPlane::interpolate
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
