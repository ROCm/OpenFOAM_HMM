/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2016-2017 Wikki Ltd
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

\*----------------------------------------------------------------------------*/

#include "volSurfaceMapping.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::volSurfaceMapping::mapToSurface
(
    const typename GeometricField<Type, fvPatchField, volMesh>::Boundary& df
) const
{
    // Grab labels for all faces in faMesh
    const labelList& faceLabels = mesh_.faceLabels();

    tmp<Field<Type>> tresult
    (
        new Field<Type>
        (
            faceLabels.size(),
            pTraits<Type>::zero
        )
    );
    Field<Type>& result = tresult.ref();

    // Get reference to volume mesh
    const polyMesh& pMesh = mesh_();
    const polyBoundaryMesh& bm = pMesh.boundaryMesh();

    label patchID, faceID;

    // Grab droplet cloud source by identifying patch and face
    forAll(faceLabels, i)
    {
        // Escape if face is beyond active faces, eg belongs to a face zone
        if (faceLabels[i] < pMesh.nFaces())
        {
            patchID = bm.whichPatch(faceLabels[i]);
            faceID = bm[patchID].whichFace(faceLabels[i]);

            result[i] = df[patchID][faceID];
        }
    }

    return tresult;
}


template<class Type>
void Foam::volSurfaceMapping::mapToVolume
(
    const GeometricField<Type, faPatchField, areaMesh>& af,
    typename GeometricField<Type, fvPatchField, volMesh>::Boundary& bf
) const
{
    // Grab labels for all faces in faMesh
    const labelList& faceLabels = mesh_.faceLabels();

    // Get reference to volume mesh
    const polyMesh& pMesh = mesh_();
    const polyBoundaryMesh& bm = pMesh.boundaryMesh();

    label patchID, faceID;

    const Field<Type>& afi = af.internalField();

    forAll(faceLabels, i)
    {
        // Escape if face is beyond active faces, eg belongs to a face zone
        if (faceLabels[i] < pMesh.nFaces())
        {
            patchID = bm.whichPatch(faceLabels[i]);
            faceID = bm[patchID].whichFace(faceLabels[i]);

            bf[patchID][faceID] = afi[i];
        }
    }
}


template<class Type>
void Foam::volSurfaceMapping::mapToVolume
(
    const tmp<GeometricField<Type, faPatchField, areaMesh>>& taf,
     typename GeometricField<Type, fvPatchField, volMesh>::Boundary& bf
) const
{
    mapToVolume(taf(), bf);

    taf.clear();
}


// ************************************************************************* //
