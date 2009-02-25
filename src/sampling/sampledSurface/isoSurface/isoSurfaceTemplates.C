/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    const scalar s0,
    const Type& p0,
    const bool hasSnap0,
    const Type& snapP0,

    const scalar s1,
    const Type& p1,
    const bool hasSnap1,
    const Type& snapP1
) const
{
    scalar d = s1-s0;

    if (mag(d) > VSMALL)
    {
        scalar s = (iso_-s0)/d;

        if (hasSnap1 && s >= 0.5 && s <= 1)
        {
            return snapP1;
        }
        else if (hasSnap0 && s >= 0.0 && s <= 0.5)
        {
            return snapP0;
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
    const scalar s0,
    const Type& p0,
    const bool hasSnap0,
    const Type& snapP0,

    const scalar s1,
    const Type& p1,
    const bool hasSnap1,
    const Type& snapP1,

    const scalar s2,
    const Type& p2,
    const bool hasSnap2,
    const Type& snapP2,

    const scalar s3,
    const Type& p3,
    const bool hasSnap3,
    const Type& snapP3,

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
            points.append
            (
                generatePoint(s0,p0,hasSnap0,snapP0,s1,p1,hasSnap1,snapP1)
            );
            points.append
            (
                generatePoint(s0,p0,hasSnap0,snapP0,s2,p2,hasSnap2,snapP2)
            );
            points.append
            (
                generatePoint(s0,p0,hasSnap0,snapP0,s3,p3,hasSnap3,snapP3)
            );
        break;

        case 0x0D:
        case 0x02:
            points.append
            (
                generatePoint(s1,p1,hasSnap1,snapP1,s0,p0,hasSnap0,snapP0)
            );
            points.append
            (
                generatePoint(s1,p1,hasSnap1,snapP1,s3,p3,hasSnap3,snapP3)
            );
            points.append
            (
                generatePoint(s1,p1,hasSnap1,snapP1,s2,p2,hasSnap2,snapP2)
            );
        break;

        case 0x0C:
        case 0x03:
        {
            Type tp1 =
                generatePoint(s0,p0,hasSnap0,snapP0,s2,p2,hasSnap2,snapP2);
            Type tp2 =
                generatePoint(s1,p1,hasSnap1,snapP1,s3,p3,hasSnap3,snapP3);

            points.append
            (
                generatePoint(s0,p0,hasSnap0,snapP0,s3,p3,hasSnap3,snapP3)
            );
            points.append(tp1);
            points.append(tp2);
            points.append(tp2);
            points.append
            (
                generatePoint(s1,p1,hasSnap1,snapP1,s2,p2,hasSnap2,snapP2)
            );
            points.append(tp1);
        }
        break;

        case 0x0B:
        case 0x04:
        {
            points.append
            (
                generatePoint(s2,p2,hasSnap2,snapP2,s0,p0,hasSnap0,snapP0)
            );
            points.append
            (
                generatePoint(s2,p2,hasSnap2,snapP2,s1,p1,hasSnap1,snapP1)
            );
            points.append
            (
                generatePoint(s2,p2,hasSnap2,snapP2,s3,p3,hasSnap3,snapP3)
            );
        }
        break;

        case 0x0A:
        case 0x05:
        {
            Type tp0 =
                generatePoint(s0,p0,hasSnap0,snapP0,s1,p1,hasSnap1,snapP1);
            Type tp1 =
                generatePoint(s2,p2,hasSnap2,snapP2,s3,p3,hasSnap3,snapP3);

            points.append(tp0);
            points.append(tp1);
            points.append
            (
                generatePoint(s0,p0,hasSnap0,snapP0,s3,p3,hasSnap3,snapP3)
            );
            points.append(tp0);
            points.append
            (
                generatePoint(s1,p1,hasSnap1,snapP1,s2,p2,hasSnap2,snapP2)
            );
            points.append(tp1);
        }
        break;

        case 0x09:
        case 0x06:
        {
            Type tp0 =
                generatePoint(s0,p0,hasSnap0,snapP0,s1,p1,hasSnap1,snapP1);
            Type tp1 =
                generatePoint(s2,p2,hasSnap2,snapP2,s3,p3,hasSnap3,snapP3);

            points.append(tp0);
            points.append
            (
                generatePoint(s1,p1,hasSnap1,snapP1,s3,p3,hasSnap3,snapP3)
            );
            points.append(tp1);
            points.append(tp0);
            points.append
            (
                generatePoint(s0,p0,hasSnap0,snapP0,s2,p2,hasSnap2,snapP2)
            );
            points.append(tp1);
        }
        break;

        case 0x07:
        case 0x08:
            points.append
            (
                generatePoint(s3,p3,hasSnap3,snapP3,s0,p0,hasSnap0,snapP0)
            );
            points.append
            (
                generatePoint(s3,p3,hasSnap3,snapP3,s2,p2,hasSnap2,snapP2)
            );
            points.append
            (
                generatePoint(s3,p3,hasSnap3,snapP3,s1,p1,hasSnap1,snapP1)
            );
        break;
    }
}


template<class Type>
Foam::label Foam::isoSurface::generateFaceTriPoints
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
    const bool hasNeiSnap,
    const Type& neiSnapPt,

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
            pVals[pointI],
            pCoords[pointI],
            snappedPoint[pointI] != -1,
            (
                snappedPoint[pointI] != -1
              ? snappedPoints[snappedPoint[pointI]]
              : pTraits<Type>::zero
            ),

            pVals[nextPointI],
            pCoords[nextPointI],
            snappedPoint[nextPointI] != -1,
            (
                snappedPoint[nextPointI] != -1
              ? snappedPoints[snappedPoint[nextPointI]]
              : pTraits<Type>::zero
            ),

            cVals[own],
            cCoords[own],
            snappedCc[own] != -1,
            (
                snappedCc[own] != -1
              ? snappedPoints[snappedCc[own]]
              : pTraits<Type>::zero
            ),

            neiVal,
            neiPt,
            hasNeiSnap,
            neiSnapPt,

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

    if
    (
        (cVals.size() != mesh_.nCells())
     || (pVals.size() != mesh_.nPoints())
     || (cCoords.size() != mesh_.nCells())
     || (pCoords.size() != mesh_.nPoints())
     || (snappedCc.size() != mesh_.nCells())
     || (snappedPoint.size() != mesh_.nPoints())
    )
    {
        FatalErrorIn("isoSurface::generateTriPoints(..)")
            << "Incorrect size." << endl
            << "mesh: nCells:" << mesh_.nCells()
            << " points:" << mesh_.nPoints() << endl
            << "cVals:" << cVals.size() << endl
            << "cCoords:" << cCoords.size() << endl
            << "snappedCc:" << snappedCc.size() << endl
            << "pVals:" << pVals.size() << endl
            << "pCoords:" << pCoords.size() << endl
            << "snappedPoint:" << snappedPoint.size() << endl
            << abort(FatalError);
    }



    // Generate triangle points

    triPoints.clear();
    triMeshCells.clear();

    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        if (faceCutType_[faceI] != NOTCUT)
        {
            generateFaceTriPoints
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
                snappedCc[nei[faceI]] != -1,
                (
                    snappedCc[nei[faceI]] != -1
                  ? snappedPoints[snappedCc[nei[faceI]]]
                  : pTraits<Type>::zero
                ),

                triPoints,
                triMeshCells
            );
        }
    }


    // Determine neighbouring snap status
    boolList neiSnapped(mesh_.nFaces()-mesh_.nInternalFaces(), false);
    List<Type> neiSnappedPoint(neiSnapped.size(), pTraits<Type>::zero);
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            label faceI = pp.start();
            forAll(pp, i)
            {
                label bFaceI = faceI-mesh_.nInternalFaces();
                label snappedIndex = snappedCc[own[faceI]];

                if (snappedIndex != -1)
                {
                    neiSnapped[bFaceI] = true;
                    neiSnappedPoint[bFaceI] = snappedPoints[snappedIndex];
                }
                faceI++;
            }
        }
    }
    syncTools::swapBoundaryFaceList(mesh_, neiSnapped, false);
    syncTools::swapBoundaryFaceList(mesh_, neiSnappedPoint, false);


    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if
        (
            isA<processorPolyPatch>(pp)
         && !refCast<const processorPolyPatch>(pp).separated()
        )
        {
            //if (refCast<const processorPolyPatch>(pp).owner())
            {
                label faceI = pp.start();

                forAll(pp, i)
                {
                    if (faceCutType_[faceI] != NOTCUT)
                    {
                        label bFaceI = faceI-mesh_.nInternalFaces();
                        if
                        (
                            neiSnapped[bFaceI]
                         && (neiSnappedPoint[bFaceI]==pTraits<Type>::zero)
                        )
                        {
                            FatalErrorIn("isoSurface::generateTriPoints(..)")
                                << "problem:" << abort(FatalError);
                        }


                        generateFaceTriPoints
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
                            neiSnapped[faceI-mesh_.nInternalFaces()],
                            neiSnappedPoint[faceI-mesh_.nInternalFaces()],

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
            // Assume zero-gradient. But what about coordinates?
            label faceI = pp.start();

            forAll(pp, i)
            {
                if (faceCutType_[faceI] != NOTCUT)
                {
                    generateFaceTriPoints
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
                        false,  // fc not snapped
                        pTraits<Type>::zero,

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
                    generateFaceTriPoints
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
                        false,  // fc not snapped
                        pTraits<Type>::zero,

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
    tmp<Field<Type> > tvalues
    (
        new Field<Type>(points().size(), pTraits<Type>::zero)
    );
    Field<Type>& values = tvalues();
    labelList nValues(values.size(), 0);

    forAll(triPoints, i)
    {
        label mergedPointI = triPointMergeMap_[i];

        if (mergedPointI >= 0)
        {
            values[mergedPointI] += triPoints[i];
            nValues[mergedPointI]++;
        }
    }

    if (debug)
    {
        Pout<< "nValues:" << values.size() << endl;
        label nMult = 0;
        forAll(nValues, i)
        {
            if (nValues[i] == 0)
            {
                FatalErrorIn("isoSurface::interpolate(..)")
                    << "point:" << i << " nValues:" << nValues[i]
                    << abort(FatalError);
            }
            else if (nValues[i] > 1)
            {
                nMult++;
            }
        }
        Pout<< "Of which mult:" << nMult << endl;
    }

    forAll(values, i)
    {
        values[i] /= scalar(nValues[i]);
    }

    return tvalues;
}


// ************************************************************************* //
