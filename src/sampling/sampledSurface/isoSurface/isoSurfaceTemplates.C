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
#include "tetMatcher.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Type Foam::isoSurface::vertexInterp
(
    const scalar iso,
    const Type& p0,
    const Type& p1,
    const scalar s0,
    const scalar s1
)
{
    scalar d = s1-s0;

    if (mag(d) > VSMALL)
    {
        return (iso-s0)/d*p1 + (s1-iso)/d*p0;
    }
    else
    {
        return 0.5*(p0+p1);
    }
}


// After "Polygonising A Scalar Field Using Tetrahedrons"
// by Paul Bourke
// Get value consistent with uncompacted triangle points.
// Given tet corner sample values s0..s3 interpolate the corresponding
// values p0..p3 to construct the surface corresponding to sample value iso.
template<class Type>
void Foam::isoSurface::vertexInterp
(
    const scalar iso,
    const scalar s0,
    const scalar s1,
    const scalar s2,
    const scalar s3,

    const Type& p0,
    const Type& p1,
    const Type& p2,
    const Type& p3,

    DynamicList<Type>& points
)
{
    int triIndex = 0;
    if (s0 < iso)
    {
        triIndex |= 1;
    }
    if (s1 < iso)
    {
        triIndex |= 2;
    }
    if (s2 < iso)
    {
        triIndex |= 4;
    }
    if (s3 < iso)
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
            points.append(vertexInterp(iso,p0,p1,s0,s1));
            points.append(vertexInterp(iso,p0,p2,s0,s2));
            points.append(vertexInterp(iso,p0,p3,s0,s3));
        break;

        case 0x0D:
        case 0x02:
            points.append(vertexInterp(iso,p1,p0,s1,s0));
            points.append(vertexInterp(iso,p1,p3,s1,s3));
            points.append(vertexInterp(iso,p1,p2,s1,s2));
        break;

        case 0x0C:
        case 0x03:
        {
            const Type tp1 = vertexInterp(iso,p0,p2,s0,s2);
            const Type tp2 = vertexInterp(iso,p1,p3,s1,s3);

            points.append(vertexInterp(iso,p0,p3,s0,s3));
            points.append(tp1);
            points.append(tp2);
            points.append(tp2);
            points.append(vertexInterp(iso,p1,p2,s1,s2));
            points.append(tp1);
        }
        break;

        case 0x0B:
        case 0x04:
        {
            points.append(vertexInterp(iso,p2,p0,s2,s0));
            points.append(vertexInterp(iso,p2,p1,s2,s1));
            points.append(vertexInterp(iso,p2,p3,s2,s3));
        }
        break;

        case 0x0A:
        case 0x05:
        {
            const Type tp0 = vertexInterp(iso,p0,p1,s0,s1);
            const Type tp1 = vertexInterp(iso,p2,p3,s2,s3);

            points.append(tp0);
            points.append(tp1);
            points.append(vertexInterp(iso,p0,p3,s0,s3));
            points.append(tp0);
            points.append(vertexInterp(iso,p1,p2,s1,s2));
            points.append(tp1);
        }
        break;

        case 0x09:
        case 0x06:
        {
            const Type tp0 = vertexInterp(iso,p0,p1,s0,s1);
            const Type tp1 = vertexInterp(iso,p2,p3,s2,s3);

            points.append(tp0);
            points.append(vertexInterp(iso,p1,p3,s1,s3));
            points.append(tp1);
            points.append(tp0);
            points.append(vertexInterp(iso,p0,p2,s0,s2));
            points.append(tp1);
        }
        break;

        case 0x07:
        case 0x08:
            points.append(vertexInterp(iso,p3,p0,s3,s0));
            points.append(vertexInterp(iso,p3,p2,s3,s2));
            points.append(vertexInterp(iso,p3,p1,s3,s1));
        break;
    }
}


template <class Type>
Foam::tmp<Foam::Field<Type> >
Foam::isoSurface::sample(const Field<Type>& vField) const
{
    return tmp<Field<Type> >(new Field<Type>(vField, meshCells()));
}


template <class Type>
Foam::tmp<Foam::Field<Type> >
Foam::isoSurface::interpolate
(
    const Field<Type>& sampleCellValues,
    const Field<Type>& samplePointValues
) const
{
    tetMatcher tet;

    DynamicList<Type> triValues;

    // Note: in same order as construction of triSurface
    label oldCellI = -1;
    forAll(meshCells_, triI)
    {
        label cellI = meshCells_[triI];

        if (cellI != oldCellI)
        {
            oldCellI = cellI;

            const cell& cFaces = mesh_.cells()[cellI];

            if (tet.isA(mesh_, cellI))
            {
                // For tets don't do cell-centre decomposition, just use the
                // tet points and values

                const face& f0 = mesh_.faces()[cFaces[0]];

                // Get the other point
                const face& f1 = mesh_.faces()[cFaces[1]];
                label oppositeI = -1;
                forAll(f1, fp)
                {
                    oppositeI = f1[fp];

                    if (findIndex(f0, oppositeI) == -1)
                    {
                        break;
                    }
                }

                vertexInterp
                (
                    iso_,
                    pointValues_[f0[0]],
                    pointValues_[f0[1]],
                    pointValues_[f0[2]],
                    pointValues_[oppositeI],

                    samplePointValues[f0[0]],
                    samplePointValues[f0[1]],
                    samplePointValues[f0[2]],
                    samplePointValues[oppositeI],

                    triValues
                );
            }
            else
            {
                forAll(cFaces, cFaceI)
                {
                    label faceI = cFaces[cFaceI];
                    const face& f = mesh_.faces()[faceI];

                    for(label fp = 1; fp < f.size() - 1; fp++)
                    {
                        vertexInterp
                        (
                            iso_,
                            pointValues_[f[0]],
                            pointValues_[f[fp]],
                            pointValues_[f[f.fcIndex(fp)]],
                            cellValues_[cellI],

                            samplePointValues[f[0]],
                            samplePointValues[f[fp]],
                            samplePointValues[f[f.fcIndex(fp)]],
                            sampleCellValues[cellI],

                            triValues
                        );
                    }
                }
            }
        }
    }

    // One value per point
    tmp<Field<Type> > tvalues(new Field<Type>(points().size()));
    Field<Type>& values = tvalues();

    forAll(triValues, i)
    {
        values[triPointMergeMap_[i]] = triValues[i];
    }

    return tvalues;
}


// ************************************************************************* //
