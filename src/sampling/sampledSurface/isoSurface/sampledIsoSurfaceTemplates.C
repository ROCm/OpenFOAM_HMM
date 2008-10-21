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
#include "isoSurface.H"
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
    const fvMesh& fvm = vField.mesh();

    if (fvm.time().timeIndex() != storedTimeIndex_)
    {
        storedTimeIndex_ = fvm.time().timeIndex();

        //- Clear any stored topo
        facesPtr_.clear();

        const volScalarField& cellFld =
            fvm.lookupObject<volScalarField>(isoField_);

        tmp<pointScalarField> pointFld
        (
            volPointInterpolation::New(fvm).interpolate(cellFld)
        );

        const isoSurface iso
        (
            fvm,
            cellFld.internalField(),
            pointFld().internalField(),
            isoVal_
        );

        const_cast<sampledIsoSurface&>(*this).triSurface::operator=(iso);
        meshCells_ = iso.meshCells();
        triPointMergeMap_ = iso.triPointMergeMap();
    }

    return tmp<Field<Type> >(new Field<Type>(vField, meshCells_));
}


template <class Type>
Foam::tmp<Foam::Field<Type> >
Foam::sampledIsoSurface::interpolateField
(
    const interpolation<Type>& interpolator
) const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    if (fvm.time().timeIndex() != storedTimeIndex_)
    {
        //- Clear any stored topo
        facesPtr_.clear();

        storedTimeIndex_ = fvm.time().timeIndex();

        const volScalarField& cellFld =
            fvm.lookupObject<volScalarField>(isoField_);

        tmp<pointScalarField> pointFld
        (
            volPointInterpolation::New(fvm).interpolate(cellFld)
        );

        const isoSurface iso
        (
            fvm,
            cellFld.internalField(),
            pointFld().internalField(),
            isoVal_
        );

        const_cast<sampledIsoSurface&>(*this).triSurface::operator=(iso);
        meshCells_ = iso.meshCells();
        triPointMergeMap_ = iso.triPointMergeMap();
    }


    // One value per point
    tmp<Field<Type> > tvalues(new Field<Type>(points().size()));
    Field<Type>& values = tvalues();

    boolList pointDone(points().size(), false);

    forAll(faces(), cutFaceI)
    {
        const face& f = faces()[cutFaceI];

        forAll(f, faceVertI)
        {
            label pointI = f[faceVertI];

            if (!pointDone[pointI])
            {
                values[pointI] = interpolator.interpolate
                (
                    points()[pointI],
                    meshCells_[cutFaceI]
                );
                pointDone[pointI] = true;
            }
        }
    }

    return tvalues;
}


// ************************************************************************* //
