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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::featureEdgeMesh::featurePoint(label ptI) const
{
    return ptI < nonFeatureStart_;
}


const Foam::vectorField& Foam::featureEdgeMesh::normals() const
{
    return normals_;
}

const Foam::vectorField& Foam::featureEdgeMesh::edgeDirections() const
{
    return edgeDirections_;
}


const Foam::labelListList& Foam::featureEdgeMesh::edgeNormals() const
{
    return edgeNormals_;
}


Foam::vectorField Foam::featureEdgeMesh::edgeNormals
(
    const labelList& edgeNormIs) const
{
    vectorField norms(edgeNormIs.size());

    forAll(edgeNormIs, i)
    {
        norms[i] = normals_[edgeNormIs[i]];
    }

    return norms;
}


Foam::vectorField Foam::featureEdgeMesh::edgeNormals(label edgeI) const
{
    return edgeNormals(edgeNormals_[edgeI]);
}


const Foam::labelListList& Foam::featureEdgeMesh::featurePointNormals() const
{
    return featurePointNormals_;
}


Foam::vectorField Foam::featureEdgeMesh::featurePointNormals(label ptI) const
{
    if (!featurePoint(ptI))
    {
        WarningIn("vectorField featureEdgeMesh::featurePointNormals")
            << "Requesting the normals of a non-feature point. "
            << "Returned zero length vectorField."
            << endl;

        return vectorField(0);
    }

    labelList featPtNormIs(featurePointNormals_[ptI]);

    vectorField norms(featPtNormIs.size());

    forAll(featPtNormIs, i)
    {
        norms[i] = normals_[featPtNormIs[i]];
    }

    return norms;
}


const Foam::labelList& Foam::featureEdgeMesh::regionEdges() const
{
    return regionEdges_;
}


Foam::featureEdgeMesh::pointStatus Foam::featureEdgeMesh::getPointStatus
(
    label ptI
) const
{
    if (ptI < concaveStart_)
    {
        return CONVEX;
    }
    else if (ptI < mixedStart_)
    {
        return CONCAVE;
    }
    else if (ptI < nonFeatureStart_)
    {
        return MIXED;
    }
    else
    {
        return NONFEATURE;
    }
}


Foam::featureEdgeMesh::edgeStatus Foam::featureEdgeMesh::getEdgeStatus
(
    label edgeI
) const
{
    if (edgeI < internalStart_)
    {
        return EXTERNAL;
    }
    else if (edgeI < flatStart_)
    {
        return INTERNAL;
    }
    else if (edgeI < openStart_)
    {
        return FLAT;
    }
    else if (edgeI < multipleStart_)
    {
        return OPEN;
    }
    else
    {
        return MULTIPLE;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //





















