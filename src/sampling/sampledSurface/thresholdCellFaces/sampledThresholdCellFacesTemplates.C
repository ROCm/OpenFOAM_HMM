/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "sampledThresholdCellFaces.H"

#include "thresholdCellFaces.H"
#include "volFieldsFwd.H"
#include "pointFields.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledThresholdCellFaces::sampleOnFaces
(
    const interpolation<Type>& sampler
) const
{
    updateGeometry(); // Recreate geometry if time has changed

    const labelList& elements = meshCells_;

    const label len = faces().size();

    auto tvalues = tmp<Field<Type>>::New(len);
    auto& values = tvalues.ref();

    const faceList& fcs = faces();
    const pointField& pts = points();

    for (label i=0; i < len; ++i)
    {
        const label celli = elements[i];
        const point pt = fcs[i].centre(pts);

        values[i] = sampler.interpolate(pt, celli);
    }

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledThresholdCellFaces::sampleOnPoints
(
    const interpolation<Type>& interpolator
) const
{
    updateGeometry(); // Recreate geometry if time has changed

    const labelList& elements = meshCells_;

    // One value per point
    auto tvalues = tmp<Field<Type>>::New(points().size(), Zero);
    auto& values = tvalues.ref();

    bitSet pointDone(points().size());

    const faceList& fcs = faces();
    const pointField& pts = points();

    forAll(fcs, i)
    {
        const face& f = fcs[i];
        const label celli = elements[i];

        for (const label pointi : f)
        {
            if (pointDone.set(pointi))
            {
                values[pointi] = interpolator.interpolate
                (
                    pts[pointi],
                    celli
                );
            }
        }
    }

    return tvalues;
}


// ************************************************************************* //
