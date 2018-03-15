/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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

#include "sampledPlane.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledPlane::sampleField
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
) const
{
    return tmp<Field<Type>>::New(vField, meshCells());
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledPlane::interpolateField
(
    const interpolation<Type>& interpolator
) const
{
    // One value per point.
    // Initialize with Zero to handle missed/degenerate faces

    tmp<Field<Type>> tvalues(new Field<Type>(points().size(), Zero));
    Field<Type>& values = tvalues.ref();

    boolList pointDone(points().size(), false);

    forAll(faces(), cutFacei)
    {
        const face& f = faces()[cutFacei];
        const label celli = meshCells()[cutFacei];

        for (const label pointi : f)
        {
            if (!pointDone[pointi])
            {
                pointDone[pointi] = true;
                values[pointi] = interpolator.interpolate
                (
                    points()[pointi],
                    celli
                );
            }
        }
    }

    return tvalues;
}


// ************************************************************************* //
