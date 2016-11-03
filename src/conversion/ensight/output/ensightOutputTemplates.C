/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

template<class Type>
void Foam::ensightOutput::writeField
(
    const char* key,
    const Field<Type>& fld,
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
                    IPstream fromSlave(Pstream::scheduled, slave);
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

                OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
                toMaster
                    << fld.component(cmpt);
            }
        }
    }
}


template<class Type>
bool Foam::ensightOutput::writePatchField
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

        const List<ensightFaces::elemType> enums =
            ensightFaces::elemEnum.enums();

        forAllConstIter(List<ensightFaces::elemType>, enums, iter)
        {
            const ensightFaces::elemType& what = *iter;

            writeField
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
bool Foam::ensightOutput::writeVolField
(
    const Field<Type>& vf,
    const ensightCells& ensCells,
    ensightFile& os,
    const bool deprecatedOrder
)
{
    if (ensCells.total())
    {
        if (Pstream::master())
        {
            os.beginPart(ensCells.index());
        }

        if (deprecatedOrder)
        {
            // element ordering used in older versions
            ensightCells::elemType oldOrder[5] =
            {
                ensightCells::HEXA8,
                ensightCells::PENTA6,
                ensightCells::PYRAMID5,
                ensightCells::TETRA4,
                ensightCells::NFACED
            };

            for (int i=0; i < 5; ++i)
            {
                const ensightCells::elemType& what = oldOrder[i];

                writeField
                (
                    ensightCells::key(what),
                    Field<Type>(vf, ensCells.cellIds(what)),
                    os
                );
            }
        }
        else
        {
            const List<ensightCells::elemType> enums =
                ensightCells::elemEnum.enums();

            forAllConstIter(List<ensightCells::elemType>, enums, iter)
            {
                const ensightCells::elemType& what = *iter;

                writeField
                (
                    ensightCells::key(what),
                    Field<Type>(vf, ensCells.cellIds(what)),
                    os
                );
            }
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
        writeVolField(vf, meshCells, os, ensMesh.deprecatedOrder());
    }

    //
    // write patches
    // use sortedToc for extra safety
    //
    const labelList patchIds = patchLookup.sortedToc();
    forAll(patchIds, listi)
    {
        const label patchId   = patchIds[listi];
        const word& patchName = patchLookup[listi];
        const ensightFaces& ensFaces = patchFaces[patchName];

        writePatchField
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

            writePatchField(values, ensFaces, os);
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

        writeField
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

        writeField
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

        writeField
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
