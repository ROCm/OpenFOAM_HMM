/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

#include "sampledIsoSurface.H"
#include "volFieldsFwd.H"
#include "pointFields.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class Type>
Foam::tmp<Foam::Field<Type> >
Foam::sampledIsoSurface::sampleField
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
) const
{
    // Recreate geometry if time has changed
    createGeometry();
    return tmp<Field<Type> >(new Field<Type>(vField, surface().meshCells()));
}


template <class Type>
Foam::tmp<Foam::Field<Type> >
Foam::sampledIsoSurface::interpolateField
(
    const interpolation<Type>& interpolator
) const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    // Get fields to sample. Assume volPointInterpolation!
    const GeometricField<Type, fvPatchField, volMesh>& volFld =
        interpolator.psi();

    tmp<GeometricField<Type, pointPatchField, pointMesh> > tpointFld
    (
        volPointInterpolation::New(fvm).interpolate(volFld)
    );

    const GeometricField<Type, pointPatchField, pointMesh>& pointFld =
        tpointFld();

    // Get pointers to sampling field (both original and interpolated one)
    getIsoFields();
    // Recreate geometry if time has changed
    createGeometry();

    // Sample.
    return surface().interpolate
    (
        *volFieldPtr_,
        *pointFieldPtr_,
        volFld,
        pointFld
    );
}


// ************************************************************************* //
