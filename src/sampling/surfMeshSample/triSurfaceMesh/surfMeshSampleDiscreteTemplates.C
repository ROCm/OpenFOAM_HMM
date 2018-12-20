/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "surfMeshSampleDiscrete.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::surfMeshSampleDiscrete::sampleOnFaces
(
    const interpolation<Type>& sampler
) const
{
    if (onBoundary())
    {
        return SurfaceSource::sampleOnFaces(sampler);
    }

    // Sample cells
    return surfMeshSample::sampleOnFaces
    (
        sampler,
        sampleElements(),  //< sampleElements == meshCells
        surface().faces(),
        surface().points()
    );
}


template<class Type>
bool Foam::surfMeshSampleDiscrete::sampleType
(
    const word& fieldName,
    const word& sampleScheme
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    const auto* volFldPtr =
        SurfaceSource::mesh().findObject<VolFieldType>(fieldName);

    if (!volFldPtr)
    {
        return false;
    }

    auto samplerPtr = interpolation<Type>::New(sampleScheme, *volFldPtr);

    getOrCreateSurfField<Type>(*volFldPtr).field() =
        sampleOnFaces(*samplerPtr);

    return true;
}


// ************************************************************************* //
