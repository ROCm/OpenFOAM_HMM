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

#include "isoSurface.H"
#include "polyMesh.H"
#include "syncTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Type Foam::isoSurface::generatePoint
(
    const DynamicList<Type>& snappedPoints,

    const scalar s0,
    const Type& p0,
    const label p0Index,

    const scalar s1,
    const Type& p1,
    const label p1Index
) const
{
    scalar d = s1-s0;

    if (mag(d) > VSMALL)
    {
        scalar s = (iso_-s0)/d;

        if (s >= 0.5 && s <= 1 && p1Index != -1)
        {
            return snappedPoints[p1Index];
        }
        else if (s >= 0.0 && s <= 0.5 && p0Index != -1)
        {
            return snappedPoints[p0Index];
        }
        else
        {
            return s*p1 + (1.0-s)*p0;
        }
    }
    else
    {
        scalar s = 0.4999;

        return s*p1 + (1.0-s)*p0;
    }
}


template<class Type>
void Foam::isoSurface::generateTriPoints
(
    const DynamicList<Type>& snapped,

    const scalar s0,
    const Type& p0,
    const label p0Index,

    const scalar s1,
    const Type& p1,
    const label p1Index,

    const scalar s2,
    const Type& p2,
    const label p2Index,

    const scalar s3,
    const Type& p3,
    const label p3Index,

    DynamicList<Type>& points
) const
{
    int triIndex = 0;
    if (s0 < iso_)
    {
        triIndex |= 1;
    }
    if (s1 < iso_)
    {
        triIndex |= 2;
    }
    if (s2 < iso_)
    {
        triIndex |= 4;
    }
    if (s3 < iso_)
    {
        triIndex |= 8;
    }

    /* Form the vertices of the triangles for each case */
    switch (triIndex)
    {
        case 0x00:
        case 0x0F:
        break;

        case 0x0E:
        case 0x01:
            points.append(generatePoint(snapped,s0,p0,p0Index,s1,p1,p1Index));
            points.append(generatePoint(snapped,s0,p0,p0Index,s2,p2,p2Index));
            points.append(generatePoint(snapped,s0,p0,p0Index,s3,p3,p3Index));
        break;

        case 0x0D:
        case 0x02:
            points.append(generatePoint(snapped,s1,p1,p1Index,s0,p0,p0Index));
            points.append(generatePoint(snapped,s1,p1,p1Index,s3,p3,p3Index));
            points.append(generatePoint(snapped,s1,p1,p1Index,s2,p2,p2Index));
        break;

        case 0x0C:
        case 0x03:
        {
            Type tp1 = generatePoint(snapped,s0,p0,p0Index,s2,p2,p2Index);
            Type tp2 = generatePoint(snapped,s1,p1,p1Index,s3,p3,p3Index);

            points.append(generatePoint(snapped,s0,p0,p0Index,s3,p3,p3Index));
            points.append(tp1);
            points.append(tp2);
            points.append(tp2);
            points.append(generatePoint(snapped,s1,p1,p1Index,s2,p2,p2Index));
            points.append(tp1);
        }
        break;

        case 0x0B:
        case 0x04:
        {
            points.append(generatePoint(snapped,s2,p2,p2Index,s0,p0,p0Index));
            points.append(generatePoint(snapped,s2,p2,p2Index,s1,p1,p1Index));
            points.append(generatePoint(snapped,s2,p2,p2Index,s3,p3,p3Index));
        }
        break;

        case 0x0A:
        case 0x05:
        {
            Type tp0 = generatePoint(snapped,s0,p0,p0Index,s1,p1,p1Index);
            Type tp1 = generatePoint(snapped,s2,p2,p2Index,s3,p3,p3Index);

            points.append(tp0);
            points.append(tp1);
            points.append(generatePoint(snapped,s0,p0,p0Index,s3,p3,p3Index));
            points.append(tp0);
            points.append(generatePoint(snapped,s1,p1,p1Index,s2,p2,p2Index));
            points.append(tp1);
        }
        break;

        case 0x09:
        case 0x06:
        {
            Type tp0 = generatePoint(snapped,s0,p0,p0Index,s1,p1,p1Index);
            Type tp1 = generatePoint(snapped,s2,p2,p2Index,s3,p3,p3Index);

            points.append(tp0);
            points.append(generatePoint(snapped,s1,p1,p1Index,s3,p3,p3Index));
            points.append(tp1);
            points.append(tp0);
            points.append(generatePoint(snapped,s0,p0,p0Index,s2,p2,p2Index));
            points.append(tp1);
        }
        break;

        case 0x07:
        case 0x08:
            points.append(generatePoint(snapped,s3,p3,p3Index,s0,p0,p0Index));
            points.append(generatePoint(snapped,s3,p3,p3Index,s2,p2,p2Index));
            points.append(generatePoint(snapped,s3,p3,p3Index,s1,p1,p1Index));
        break;
    }
}


template<class Type>
Foam::label Foam::isoSurface::generateTriPoints
(
    const volScalarField& cVals,
    const scalarField& pVals,

    const GeometricField<Type, fvPatchField, volMesh>& cCoords,
    const Field<Type>& pCoords,

    const DynamicList<Type>& snappedPoints,
    const labelList& snappedCc,
    const labelList& snappedPoint,
    const label faceI,

    const scalar neiVal,
    const Type& neiPt,
    const label neiSnap,

    DynamicList<Type>& triPoints,
    DynamicList<label>& triMeshCells
) const
{
    label own = mesh_.faceOwner()[faceI];

    label oldNPoints = triPoints.size();

    const face& f = mesh_.faces()[faceI];

    forAll(f, fp)
    {
        label pointI = f[fp];
        label nextPointI = f[f.fcIndex(fp)];

        generateTriPoints
        (
            snappedPoints,

            pVals[pointI],
            pCoords[pointI],
            snappedPoint[pointI],

            pVals[nextPointI],
            pCoords[nextPointI],
            snappedPoint[nextPointI],

            cVals[own],
            cCoords[own],
            snappedCc[own],

            neiVal,
            neiPt,
            neiSnap,

            triPoints
        );
    }

    // Every three triPoints is a triangle
    label nTris = (triPoints.size()-oldNPoints)/3;
    for (label i = 0; i < nTris; i++)
    {
        triMeshCells.append(own);
    }

    return nTris;
}


template<class Type>
void Foam::isoSurface::generateTriPoints
(
    const volScalarField& cVals,
    const scalarField& pVals,

    const GeometricField<Type, fvPatchField, volMesh>& cCoords,
    const Field<Type>& pCoords,

    const DynamicList<Type>& snappedPoints,
    const labelList& snappedCc,
    const labelList& snappedPoint,

    DynamicList<Type>& triPoints,
    DynamicList<label>& triMeshCells
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();


    // Determine neighbouring snap status
    labelList neiSnappedCc(mesh_.nFaces()-mesh_.nInternalFaces(), -1);
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            label faceI = pp.start();
            forAll(pp, i)
            {
                neiSnappedCc[faceI-mesh_.nInternalFaces()] =
                    snappedCc[own[faceI]];
                faceI++;
            }
        }
    }
    syncTools::swapBoundaryFaceList(mesh_, neiSnappedCc, false);



    // Generate triangle points

    triPoints.clear();
    triMeshCells.clear();

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        if (faceCutType_[faceI] != NOTCUT)
        {
            generateTriPoints
            (
                cVals,
                pVals,

                cCoords,
                pCoords,

                snappedPoints,
                snappedCc,
                snappedPoint,
                faceI,

                cVals[nei[faceI]],
                cCoords[nei[faceI]],
                snappedCc[nei[faceI]],

                triPoints,
                triMeshCells
            );
        }
    }


    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            if (refCast<const processorPolyPatch>(pp).owner())
            {
                label faceI = pp.start();

                forAll(pp, i)
                {
                    if (faceCutType_[faceI] != NOTCUT)
                    {
                        generateTriPoints
                        (
                            cVals,
                            pVals,

                            cCoords,
                            pCoords,

                            snappedPoints,
                            snappedCc,
                            snappedPoint,
                            faceI,

                            cVals.boundaryField()[patchI][i],
                            cCoords.boundaryField()[patchI][i],
                            neiSnappedCc[faceI-mesh_.nInternalFaces()],

                            triPoints,
                            triMeshCells
                        );
                    }
                    faceI++;
                }
            }
        }
        else if (isA<emptyPolyPatch>(pp))
        {
            // Assume zero-gradient.
            label faceI = pp.start();

            forAll(pp, i)
            {
                if (faceCutType_[faceI] != NOTCUT)
                {
                    generateTriPoints
                    (
                        cVals,
                        pVals,

                        cCoords,
                        pCoords,

                        snappedPoints,
                        snappedCc,
                        snappedPoint,
                        faceI,

                        cVals[own[faceI]],
                        cCoords.boundaryField()[patchI][i],
                        -1, // fc not snapped

                        triPoints,
                        triMeshCells
                    );
                }
                faceI++;
            }
        }
        else
        {
            label faceI = pp.start();

            forAll(pp, i)
            {
                if (faceCutType_[faceI] != NOTCUT)
                {
                    generateTriPoints
                    (
                        cVals,
                        pVals,

                        cCoords,
                        pCoords,

                        snappedPoints,
                        snappedCc,
                        snappedPoint,
                        faceI,

                        cVals.boundaryField()[patchI][i],
                        cCoords.boundaryField()[patchI][i],
                        -1, // fc not snapped

                        triPoints,
                        triMeshCells
                    );
                }
                faceI++;
            }
        }
    }

    triPoints.shrink();
    triMeshCells.shrink();
}


//template <class Type>
//Foam::tmp<Foam::Field<Type> >
//Foam::isoSurface::sample(const Field<Type>& vField) const
//{
//    return tmp<Field<Type> >(new Field<Type>(vField, meshCells()));
//}
//
//
template <class Type>
Foam::tmp<Foam::Field<Type> >
Foam::isoSurface::interpolate
(
    const volScalarField& cVals,
    const scalarField& pVals,
    const GeometricField<Type, fvPatchField, volMesh>& cCoords,
    const Field<Type>& pCoords
) const
{
    DynamicList<Type> triPoints(nCutCells_);
    DynamicList<label> triMeshCells(nCutCells_);

    // Dummy snap data
    DynamicList<Type> snappedPoints;
    labelList snappedCc(mesh_.nCells(), -1);
    labelList snappedPoint(mesh_.nPoints(), -1);

    generateTriPoints
    (
        cVals,
        pVals,

        cCoords,
        pCoords,

        snappedPoints,
        snappedCc,
        snappedPoint,

        triPoints,
        triMeshCells
    );


    // One value per point
    tmp<Field<Type> > tvalues(new Field<Type>(points().size()));
    Field<Type>& values = tvalues();

    forAll(triPoints, i)
    {
        label mergedPointI = triPointMergeMap_[i];

        if (mergedPointI >= 0)
        {
            values[mergedPointI] = triPoints[i];
        }
    }

    return tvalues;
}


// ************************************************************************* //
