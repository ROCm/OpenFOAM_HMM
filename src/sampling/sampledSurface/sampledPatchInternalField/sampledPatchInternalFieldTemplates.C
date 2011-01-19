/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "sampledPatchInternalField.H"
#include "interpolationCellPoint.H"
#include "PrimitivePatchInterpolation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class Type>
Foam::tmp<Foam::Field<Type> >
Foam::sampledPatchInternalField::sampleField
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
) const
{
    const mapDistribute& distMap = map();

    // One value per face
    tmp<Field<Type> > tvalues(new Field<Type>(patchFaceLabels().size()));
    Field<Type>& values = tvalues();

    if (patchIndex() != -1)
    {
        Field<Type> interpVals = vField.internalField();
        distMap.distribute(interpVals);

        forAll(patchFaceLabels(), elemI)
        {
            values[elemI] = interpVals[patchFaceLabels()[elemI]];
        }
    }

    return tvalues;
}


template <class Type>
Foam::tmp<Foam::Field<Type> >
Foam::sampledPatchInternalField::interpolateField
(
    const interpolation<Type>& interpolator
) const
{
    // One value per vertex

    if (patchIndex() != -1)
    {
        // See directMappedFixedValueFvPatchField
        const mapDistribute& distMap = map();

        const polyPatch& pp = mesh().boundaryMesh()[patchIndex()];

        // Send back sample points to processor that holds the cell.
        // Mark cells with point::max so we know which ones we need
        // to interpolate (since expensive).
        vectorField samples(samplePoints());
        distMap.reverseDistribute(mesh().nCells(), point::max, samples);

        Field<Type> patchVals(mesh().nCells());

        forAll(samples, cellI)
        {
            if (samples[cellI] != point::max)
            {
                patchVals[cellI] = interpolator.interpolate
                (
                    samples[cellI],
                    cellI
                );
            }
        }

        distMap.distribute(patchVals);

        // Now patchVals holds the interpolated data in patch face order.
        // Interpolate to points. Note: points are patch.localPoints() so
        // can use standard interpolation

        return PrimitivePatchInterpolation<primitivePatch>
        (
           pp
        ).faceToPointInterpolate(patchVals);
    }
    else
    {
        return tmp<Field<Type> >(new Field<Type>(points().size()));
    }
}


// ************************************************************************* //
