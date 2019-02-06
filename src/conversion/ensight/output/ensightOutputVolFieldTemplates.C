/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "ensightOutputVolField.H"
#include "ensightMesh.H"

#include "fvMesh.H"
#include "globalIndex.H"
#include "volPointInterpolation.H"
#include "interpolation.H"
#include "linear.H"
#include "processorFvPatch.H"
#include "uindirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * *  Detail * * * * * * * * * * * * * * * * * //

template<class Type>
bool Foam::ensightOutput::Detail::writeVolField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const ensightMesh& ensMesh,
    ensightFile& os
)
{
    constexpr bool parallel = true;

    const fvMesh& mesh = ensMesh.mesh();
    const ensightCells& meshCells = ensMesh.meshCells();
    const Map<word>&  patchLookup = ensMesh.patches();
    const HashTable<ensightFaces>& patchFaces = ensMesh.boundaryPatchFaces();
    const HashTable<ensightFaces>&  zoneFaces = ensMesh.faceZoneFaces();

    //
    // write internalMesh, unless patch-selection was requested
    //
    if (ensMesh.useInternalMesh())
    {
        Detail::writeCellField(vf, meshCells, os, parallel);
    }

    //
    // write patches
    // use sortedToc for extra safety
    //
    const labelList patchIds = patchLookup.sortedToc();
    for (const label patchId : patchIds)
    {
        const word& patchName = patchLookup[patchId];
        const ensightFaces& part = patchFaces[patchName];

        writeFaceField
        (
            vf.boundaryField()[patchId],
            part,
            os,
            parallel
        );
    }


    //
    // write faceZones, if requested
    // use sortedToc for extra safety
    //
    const wordList zoneNames = zoneFaces.sortedToc();
    if (!zoneNames.empty())
    {
        // Interpolates cell values to faces - needed only when exporting
        // faceZones...
        GeometricField<Type, fvsPatchField, surfaceMesh> sf
        (
            Foam::linearInterpolate(vf)
        );

        // flat boundary field
        // similar to volPointInterpolation::flatBoundaryField()

        Field<Type> flat(mesh.nBoundaryFaces(), Zero);

        const fvBoundaryMesh& bm = mesh.boundary();
        forAll(vf.boundaryField(), patchi)
        {
            const polyPatch& pp = bm[patchi].patch();
            const auto& bf = vf.boundaryField()[patchi];

            if (isA<processorFvPatch>(bm[patchi]))
            {
                // Use average value for processor faces
                // own cell value = patchInternalField
                // nei cell value = evaluated boundary values
                SubList<Type>
                (
                    flat,
                    bf.size(),
                    pp.offset()
                ) = (0.5 * (bf.patchInternalField() + bf));
            }
            else if (!isA<emptyFvPatch>(bm[patchi]))
            {
                SubList<Type>
                (
                    flat,
                    bf.size(),
                    pp.offset()
                ) = bf;
            }
        }

        for (const word& zoneName : zoneNames)
        {
            const ensightFaces& part = zoneFaces[zoneName];

            // Field (local size)
            Field<Type> values(part.size());

            // Loop over face ids to store the needed field values
            // - internal faces use linear interpolation
            // - boundary faces use the corresponding patch value
            forAll(part, i)
            {
                const label faceId = part[i];
                values[i] =
                (
                    mesh.isInternalFace(faceId)
                  ? sf[faceId]
                  : flat[faceId - mesh.nInternalFaces()]
                );
            }

            // The field is already copied in the proper order
            // - just need its corresponding sub-fields
            Detail::writeFaceSubField(values, part, os, parallel);
        }
    }

    return true;
}


template<class Type>
bool Foam::ensightOutput::Detail::writePointField
(
    const GeometricField<Type, pointPatchField, pointMesh>& pf,
    const ensightMesh& ensMesh,
    ensightFile& os
)
{
    constexpr bool parallel = true;

    const fvMesh& mesh = ensMesh.mesh();
    const Map<word>& patchLookup  = ensMesh.patches();

    const HashTable<ensightFaces>& patchFaces = ensMesh.boundaryPatchFaces();
    const HashTable<ensightFaces>&  zoneFaces = ensMesh.faceZoneFaces();

    //
    // write internalMesh, unless patch-selection was requested
    //
    if (ensMesh.useInternalMesh())
    {
        if (Pstream::master())
        {
            os.beginPart(0); // 0 = internalMesh
        }

        Detail::writeFieldComponents
        (
            "coordinates",
            Field<Type>(pf.internalField(), ensMesh.uniquePointMap()),
            os,
            parallel
        );
    }

    //
    // write patches
    // use sortedToc for extra safety
    //
    const labelList patchIds = patchLookup.sortedToc();
    for (const label patchId : patchIds)
    {
        const word& patchName = patchLookup[patchId];
        const ensightFaces& part = patchFaces[patchName];

        const fvPatch& p = mesh.boundary()[patchId];

        // Renumber the patch points/faces into unique points
        labelList pointToGlobal;
        labelList uniqueMeshPointLabels;
        autoPtr<globalIndex> globalPointsPtr =
            mesh.globalData().mergePoints
            (
                p.patch().meshPoints(),
                p.patch().meshPointMap(),
                pointToGlobal,
                uniqueMeshPointLabels
            );

        if (Pstream::master())
        {
            os.beginPart(part.index());
        }

        Detail::writeFieldComponents
        (
            "coordinates",
            Field<Type>(pf.internalField(), uniqueMeshPointLabels),
            os,
            parallel
        );
    }

    //
    // write faceZones, if requested
    //
    const wordList zoneNames = zoneFaces.sortedToc();
    for (const word& zoneName : zoneNames)
    {
        const ensightFaces& part = zoneFaces[zoneName];

        uindirectPrimitivePatch p
        (
            UIndirectList<face>
            (
                mesh.faces(),
                part.faceIds()
            ),
            mesh.points()
        );

        // Renumber the patch points/faces into unique points
        labelList pointToGlobal;
        labelList uniqueMeshPointLabels;
        autoPtr<globalIndex> globalPointsPtr =
            mesh.globalData().mergePoints
            (
                p.meshPoints(),
                p.meshPointMap(),
                pointToGlobal,
                uniqueMeshPointLabels
            );

        if (Pstream::master())
        {
            os.beginPart(part.index());
        }

        Detail::writeFieldComponents
        (
            "coordinates",
            Field<Type>(pf.internalField(), uniqueMeshPointLabels),
            os,
            parallel
        );
    }

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
bool Foam::ensightOutput::writeVolField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const ensightMesh& ensMesh,
    ensightFile& os,
    const bool nodeValues
)
{
    if (nodeValues)
    {
        tmp<GeometricField<Type, pointPatchField, pointMesh>> pfld
        (
            volPointInterpolation::New(vf.mesh()).interpolate(vf)
        );
        pfld.ref().checkOut();
        pfld.ref().rename(vf.name());

        return Detail::writePointField<Type>(pfld, ensMesh, os);
    }

    return Detail::writeVolField<Type>(vf, ensMesh, os);
}


// * * * * * * * * * * * * * * * *  Serial * * * * * * * * * * * * * * * * * //

template<class Type>
bool Foam::ensightOutput::Serial::writeVolField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const ensightPartFaces& part,
    ensightFile& os
)
{
    const label patchi = part.patchIndex();

    if (patchi >= 0 && patchi < vf.boundaryField().size())
    {
        return ensightOutput::Detail::writeFaceField
        (
            vf.boundaryField()[patchi],
            part,
            os,
            false // serial
        );
    }

    return false;
}


template<class Type>
bool Foam::ensightOutput::Serial::writeVolField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const ensightPartCells& part,
    ensightFile& os
)
{
    return ensightOutput::Detail::writeCellField
    (
        vf.internalField(),
        part,
        os,
        false // serial
    );
}


template<class Type>
bool Foam::ensightOutput::Serial::writeVolField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const ensightParts& list,
    ensightFile& os
)
{
    for (const ensightPart& part : list)
    {
        const ensightPartFaces* fptr = isA<ensightPartFaces>(part);

        if (fptr)
        {
            Serial::writeVolField(vf, *fptr, os);
            continue;
        }

        const ensightPartCells* cptr = isA<ensightPartCells>(part);

        if (cptr)
        {
            Serial::writeVolField(vf, *cptr, os);
            continue;
        }
    }

    return true;
}


// ************************************************************************* //
