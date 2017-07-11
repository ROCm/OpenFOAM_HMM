/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "ensightFile.H"
#include "ensightOutput.H"
#include "ensightPTraits.H"

#include "fvMesh.H"
#include "volFields.H"
#include "IOField.H"
#include "OFstream.H"
#include "IOmanip.H"
#include "Time.H"
#include "volPointInterpolation.H"
#include "globalIndex.H"
#include "uindirectPrimitivePatch.H"
#include "interpolation.H"
#include "linear.H"

// * * * * * * * * * * Static Private Member Functions * * * * * * * * * * * //

template<template<typename> class FieldContainer, class Type>
void Foam::ensightOutput::writeFieldContent
(
    const char* key,
    const FieldContainer<Type>& fld,
    ensightFile& os
)
{
    if (returnReduce(fld.size(), sumOp<label>()) > 0)
    {
        if (Pstream::master())
        {
            os.writeKeyword(key);

            for (direction d=0; d < pTraits<Type>::nComponents; ++d)
            {
                const label cmpt = ensightPTraits<Type>::componentOrder[d];

                os.writeList(fld.component(cmpt));

                for (int slave=1; slave<Pstream::nProcs(); ++slave)
                {
                    IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
                    scalarField received(fromSlave);
                    os.writeList(received);
                }
            }
        }
        else
        {
            for (direction d=0; d < pTraits<Type>::nComponents; ++d)
            {
                const label cmpt = ensightPTraits<Type>::componentOrder[d];

                OPstream toMaster
                (
                    Pstream::commsTypes::scheduled,
                    Pstream::masterNo()
                );

                toMaster
                    << fld.component(cmpt);
            }
        }
    }
}


template<class Type>
bool Foam::ensightOutput::writeFaceField
(
    const Field<Type>& pf,
    const ensightFaces& ensFaces,
    Foam::ensightFile& os
)
{
    if (ensFaces.total())
    {
        if (Pstream::master())
        {
            os.beginPart(ensFaces.index());
        }

        for (label typei=0; typei < ensightFaces::nTypes; ++typei)
        {
            const ensightFaces::elemType what =
                ensightFaces::elemType(typei);

            writeFieldContent
            (
                ensightFaces::key(what),
                Field<Type>(pf, ensFaces.faceIds(what)),
                os
            );
        }

        return true;
    }
    else
    {
        return false;
    }
}


template<class Type>
bool Foam::ensightOutput::writeFaceSubField
(
    const Field<Type>& pf,
    const ensightFaces& ensFaces,
    Foam::ensightFile& os
)
{
    if (ensFaces.total())
    {
        if (Pstream::master())
        {
            os.beginPart(ensFaces.index());
        }

        label start = 0; // start of sublist
        for (label typei=0; typei < ensightFaces::nTypes; ++typei)
        {
            const ensightFaces::elemType what = ensightFaces::elemType(typei);
            const label size = ensFaces.faceIds(what).size();

            writeFieldContent
            (
                ensightFaces::key(what),
                SubField<Type>(pf, size, start),
                os
            );

            start += size; // start of next sublist
        }

        return true;
    }
    else
    {
        return false;
    }
}


template<class Type>
bool Foam::ensightOutput::writeCellField
(
    const Field<Type>& vf,
    const ensightCells& ensCells,
    ensightFile& os
)
{
    if (ensCells.total())
    {
        if (Pstream::master())
        {
            os.beginPart(ensCells.index());
        }

        for (label typei=0; typei < ensightCells::nTypes; ++typei)
        {
            const ensightCells::elemType what = ensightCells::elemType(typei);

            writeFieldContent
            (
                ensightCells::key(what),
                Field<Type>(vf, ensCells.cellIds(what)),
                os
            );
        }

        return true;
    }
    else
    {
        return false;
    }
}


template<class Type>
bool Foam::ensightOutput::writeField
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const ensightMesh& ensMesh,
    ensightFile& os
)
{
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
        writeCellField(vf, meshCells, os);
    }

    //
    // write patches
    // use sortedToc for extra safety
    //
    const labelList patchIds = patchLookup.sortedToc();
    forAll(patchIds, listi)
    {
        const label patchId   = patchIds[listi];
        const word& patchName = patchLookup[patchId];
        const ensightFaces& ensFaces = patchFaces[patchName];

        writeFaceField
        (
            vf.boundaryField()[patchId],
            ensFaces,
            os
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
        // as per volPointInterpolation::flatBoundaryField()

        Field<Type> flat(mesh.nFaces() - mesh.nInternalFaces());

        const fvBoundaryMesh& bm = mesh.boundary();
        forAll(vf.boundaryField(), patchI)
        {
            const polyPatch& pp = bm[patchI].patch();
            const label bFaceI = pp.start() - mesh.nInternalFaces();

            if
            (
                isA<emptyFvPatch>(bm[patchI])
             || vf.boundaryField()[patchI].coupled()
            )
            {
                SubList<Type>
                (
                    flat,
                    pp.size(),
                    bFaceI
                ) = Zero;
            }
            else
            {
                SubList<Type>
                (
                    flat,
                    vf.boundaryField()[patchI].size(),
                    bFaceI
                ) = vf.boundaryField()[patchI];
            }
        }

        forAll(zoneNames, zonei)
        {
            const word& zoneName = zoneNames[zonei];
            const ensightFaces& ensFaces = zoneFaces[zoneName];

            // field (local size)
            Field<Type> values(ensFaces.size());

            // Loop over face ids to store the needed field values
            // - internal faces use linear interpolation
            // - boundary faces use the corresponding patch value
            forAll(ensFaces, i)
            {
                label faceId = ensFaces[i];
                values[i] =
                (
                    mesh.isInternalFace(faceId)
                  ? sf[faceId]
                  : flat[faceId - mesh.nInternalFaces()]
                );
            }

            // The field is already copied in the proper order
            // - just need its corresponding sub-fields
            writeFaceSubField(values, ensFaces, os);
        }
    }

    return true;
}


template<class Type>
bool Foam::ensightOutput::ensightPointField
(
    const GeometricField<Type, pointPatchField, pointMesh>& pf,
    const ensightMesh& ensMesh,
    ensightFile& os
)
{
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

        writeFieldContent
        (
            "coordinates",
            Field<Type>(pf.internalField(), ensMesh.uniquePointMap()),
            os
        );
    }

    //
    // write patches
    // use sortedToc for extra safety
    //
    const labelList patchIds = patchLookup.sortedToc();
    forAll(patchIds, listi)
    {
        const label patchId   = patchIds[listi];
        const word& patchName = patchLookup[patchId];
        const ensightFaces& ensFaces = patchFaces[patchName];

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
            os.beginPart(ensFaces.index());
        }

        writeFieldContent
        (
            "coordinates",
            Field<Type>(pf.internalField(), uniqueMeshPointLabels),
            os
        );
    }

    //
    // write faceZones, if requested
    //
    const wordList zoneNames = zoneFaces.sortedToc();
    forAll(zoneNames, zonei)
    {
        const word& zoneName = zoneNames[zonei];
        const ensightFaces& ensFaces = zoneFaces[zoneName];

        uindirectPrimitivePatch p
        (
            UIndirectList<face>
            (
                mesh.faces(),
                ensFaces.faceIds()
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
            os.beginPart(ensFaces.index());
        }

        writeFieldContent
        (
            "coordinates",
            Field<Type>(pf.internalField(), uniqueMeshPointLabels),
            os
        );
    }

    return true;
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
bool Foam::ensightOutput::writeField
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

        return ensightPointField<Type>(pfld, ensMesh, os);
    }
    else
    {
        return writeField<Type>(vf, ensMesh, os);
    }
}


// ************************************************************************* //
