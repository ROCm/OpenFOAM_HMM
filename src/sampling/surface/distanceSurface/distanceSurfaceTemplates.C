/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 OpenCFD Ltd.
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

#include "distanceSurface.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::distanceSurface::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& cellValues,
    const Field<Type>& pointValues
) const
{
    if (isoSurfCellPtr_)
    {
        return isoSurfCellPtr_->interpolate(cellValues, pointValues);
    }
    else if (isoSurfTopoPtr_)
    {
        return isoSurfTopoPtr_->interpolate(cellValues, pointValues);
    }
    else
    {
        return isoSurfPtr_->interpolate(cellValues, pointValues);
    }
}


// ************************************************************************* //
