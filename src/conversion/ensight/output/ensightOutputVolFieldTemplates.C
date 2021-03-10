/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
#include "linear.H"
#include "volPointInterpolation.H"
#include "interpolation.H"
#include "processorFvPatch.H"
#include "DynamicField.H"

// * * * * * * * * * * * * * * * *  Detail * * * * * * * * * * * * * * * * * //

template<class Type>
bool Foam::ensightOutput::writeVolField
(
    ensightFile& os,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const ensightMesh& ensMesh
)
{
    bool parallel = Pstream::parRun();

    const fvMesh& mesh = vf.mesh();
    const polyBoundaryMesh& bmesh = mesh.boundaryMesh();

    const Map<ensightCells>& cellZoneParts = ensMesh.cellZoneParts();
    const Map<ensightFaces>& faceZoneParts = ensMesh.faceZoneParts();
    const Map<ensightFaces>& boundaryParts = ensMesh.boundaryParts();


    // Write internalMesh and cellZones - sorted by index
    for (const label zoneId : cellZoneParts.sortedToc())
    {
        const ensightCells& part = cellZoneParts[zoneId];

        ensightOutput::writeField(os, vf, part, parallel);
    }


    // Write patches - sorted by index
    for (const label patchId : boundaryParts.sortedToc())
    {
        const ensightFaces& part = boundaryParts[patchId];

        if (patchId < 0 || patchId >= bmesh.size())
        {
            // Future handling of combined patches?
            continue;
        }

        const label patchStart = bmesh[patchId].start();

        // Either use a flat boundary field for all patches,
        // or patch-local face ids

        // Operate on a copy
        ensightFaces localPart(part);

        // Change from global faceIds to patch-local
        localPart.decrFaceIds(patchStart);

        ensightOutput::writeField
        (
            os,
            vf.boundaryField()[patchId],
            localPart,
            parallel
        );
    }


    // No face zones data
    if (faceZoneParts.empty())
    {
        return true;
    }


    // Flat boundary field
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


    // Interpolate cell values to faces
    // - only needed when exporting faceZones...

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsfld
        = Foam::linearInterpolate(vf);

    const auto& sfld = tsfld();


    // Local output buffer
    label maxLen = 0;

    forAllConstIters(faceZoneParts, iter)
    {
        maxLen = max(maxLen, iter.val().size());
    }

    DynamicField<Type> values(maxLen);


    //
    // Write requested faceZones - sorted by index
    //
    for (const label zoneId : faceZoneParts.sortedToc())
    {
        const ensightFaces& part = faceZoneParts[zoneId];

        // Loop over face ids to store the needed field values
        // - internal faces use linear interpolation
        // - boundary faces use the corresponding patch value

        // Local copy of the field
        values.resize(part.size());
        values = Zero;

        auto valIter = values.begin();

        for (const label faceId : part.faceIds())
        {
            *valIter =
            (
                mesh.isInternalFace(faceId)
              ? sfld[faceId]
              : flat[faceId - mesh.nInternalFaces()]
            );

            ++valIter;
        }

        // The field is already in the proper element order
        // - just need its corresponding sub-fields
        ensightOutput::Detail::writeFaceSubField(os, values, part, parallel);
    }

    return true;
}


template<class Type>
bool Foam::ensightOutput::writePointField
(
    ensightFile& os,
    const GeometricField<Type, pointPatchField, pointMesh>& pf,
    const ensightMesh& ensMesh
)
{
    bool parallel = Pstream::parRun();

    const polyMesh& mesh = ensMesh.mesh();

    const Map<ensightCells>& cellZoneParts = ensMesh.cellZoneParts();
    const Map<ensightFaces>& faceZoneParts = ensMesh.faceZoneParts();
    const Map<ensightFaces>& boundaryParts = ensMesh.boundaryParts();

    //
    // Write internalMesh and cellZones - sorted by index
    //
    for (const label zoneId : cellZoneParts.sortedToc())
    {
        const ensightCells& part = cellZoneParts[zoneId];

        if (Pstream::master())
        {
            os.beginPart(part.index());
        }

        labelList uniquePointLabels;
        part.uniqueMeshPoints(mesh, uniquePointLabels, parallel);

        ensightOutput::Detail::writeFieldComponents
        (
            os,
            ensightFile::coordinates,
            UIndirectList<Type>(pf.internalField(), uniquePointLabels),
            parallel
        );
    }


    //
    // Write patches - sorted by index
    //
    for (const label patchId : boundaryParts.sortedToc())
    {
        const ensightFaces& part = boundaryParts[patchId];

        if (Pstream::master())
        {
            os.beginPart(part.index());
        }

        labelList uniquePointLabels;
        part.uniqueMeshPoints(mesh, uniquePointLabels, parallel);

        const auto& bfld = pf.boundaryField()[patchId];

        // Only valuePointPatchField is actually derived from Field
        const auto* vpp = isA<Field<Type>>(bfld);
        if (vpp)
        {
            // Require patch local indices, not mesh point labels.
            // But need to use polyPatch meshPointMap() to recover the
            // local indices since the ensight output will have jumbled
            // the face output order

            const polyPatch& pp = mesh.boundaryMesh()[patchId];

            for (label& pointi : uniquePointLabels)
            {
                pointi = pp.meshPointMap()[pointi];
            }

            ensightOutput::Detail::writeFieldComponents
            (
                os,
                ensightFile::coordinates,
                UIndirectList<Type>(*vpp, uniquePointLabels),
                parallel
            );
        }
        else
        {
            ensightOutput::Detail::writeFieldComponents
            (
                os,
                ensightFile::coordinates,
                UIndirectList<Type>(pf.internalField(), uniquePointLabels),
                parallel
            );
        }
    }

    //
    // Write requested faceZones - sorted by index
    //
    for (const label zoneId : faceZoneParts.sortedToc())
    {
        const ensightFaces& part = faceZoneParts[zoneId];

        if (Pstream::master())
        {
            os.beginPart(part.index());
        }

        // CAVEAT - does not properly handle valuePointPatchField,
        // uses internalField only

        {
            labelList uniquePointLabels;
            part.uniqueMeshPoints(mesh, uniquePointLabels, parallel);

            ensightOutput::Detail::writeFieldComponents
            (
                os,
                ensightFile::coordinates,
                UIndirectList<Type>(pf.internalField(), uniquePointLabels),
                parallel
            );
        }
    }

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
bool Foam::ensightOutput::writeVolField
(
    ensightFile& os,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const ensightMesh& ensMesh,
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

        return ensightOutput::writePointField<Type>(os, pfld, ensMesh);
    }

    return ensightOutput::writeVolField<Type>(os, vf, ensMesh);
}


// ************************************************************************* //
