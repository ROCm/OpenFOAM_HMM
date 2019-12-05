/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "sampledCuttingPlane.H"
#include "volFieldsFwd.H"
#include "pointFields.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledCuttingPlane::sampleOnFaces
(
    const interpolation<Type>& sampler
) const
{
    return sampledSurface::sampleOnFaces
    (
        sampler,
        meshCells(),
        faces(),
        points()
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledCuttingPlane::sampleOnPoints
(
    const interpolation<Type>& interpolator
) const
{
    // Assume volPointInterpolation for the point field!
    const auto& volFld = interpolator.psi();

    if (subMeshPtr_.valid())
    {
        auto tvolSubFld = subMeshPtr_().interpolate(volFld);
        const auto& volSubFld = tvolSubFld();

        auto tpointFld =
            volPointInterpolation::New(volSubFld.mesh()).interpolate(volSubFld);

        return this->isoSurfaceInterpolate
        (
            (average_ ? pointAverage(tpointFld())() : volSubFld),
            tpointFld()
        );
    }


    auto tpointFld =
        volPointInterpolation::New(volFld.mesh()).interpolate(volFld);

    return this->isoSurfaceInterpolate
    (
        (average_ ? pointAverage(tpointFld())() : volFld),
        tpointFld()
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledCuttingPlane::isoSurfaceInterpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& cCoords,
    const Field<Type>& pCoords
) const
{
    if (isoSurfCellPtr_)
    {
        return isoSurfCellPtr_->interpolate(cCoords, pCoords);
    }
    else if (isoSurfTopoPtr_)
    {
        return isoSurfTopoPtr_->interpolate(cCoords, pCoords);
    }
    else
    {
        return isoSurfPtr_->interpolate(cCoords, pCoords);
    }
}


// ************************************************************************* //
