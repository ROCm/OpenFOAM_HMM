/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::isoSurfaceTopo::interpolateTemplate
(
    const Field<Type>& cellCoords,
    const Field<Type>& pointCoords
) const
{
    auto tfld = tmp<Field<Type>>::New(pointToVerts_.size());
    auto& fld = tfld.ref();

    forAll(pointToVerts_, i)
    {
        scalar s0;
        Type p0;
        {
            label idx = pointToVerts_[i].first();
            if (idx < mesh_.nPoints())
            {
                // Point index
                s0 = pVals_[idx];
                p0 = pointCoords[idx];
            }
            else
            {
                // Cell index
                idx -= mesh_.nPoints();
                s0 = cVals_[idx];
                p0 = cellCoords[idx];
            }
        }

        scalar s1;
        Type p1;
        {
            label idx = pointToVerts_[i].second();
            if (idx < mesh_.nPoints())
            {
                // Point index
                s1 = pVals_[idx];
                p1 = pointCoords[idx];
            }
            else
            {
                // Cell index
                idx -= mesh_.nPoints();
                s1 = cVals_[idx];
                p1 = cellCoords[idx];
            }
        }

        const scalar d = s1-s0;
        if (mag(d) > VSMALL)
        {
            const scalar s = (iso_-s0)/d;
            fld[i] = s*p1+(1.0-s)*p0;
        }
        else
        {
            fld[i] = 0.5*(p0+p1);
        }
    }

    return tfld;
}


// ************************************************************************* //
