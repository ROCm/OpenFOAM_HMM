/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::primitiveMesh::nInternalPoints() const
{
    return nInternalPoints_;
}


inline Foam::label Foam::primitiveMesh::nPoints() const
{
    return nPoints_;
}


inline Foam::label Foam::primitiveMesh::nInternal0Edges() const
{
    // Force edge calculation
    (void)nEdges();
    return nInternal0Edges_;
}


inline Foam::label Foam::primitiveMesh::nInternal1Edges() const
{
    // Force edge calculation
    (void)nEdges();
    return nInternal1Edges_;
}


inline Foam::label Foam::primitiveMesh::nInternalEdges() const
{
    // Force edge calculation
    (void)nEdges();
    return nInternalEdges_;
}


inline Foam::label Foam::primitiveMesh::nEdges() const
{
    if (nEdges_ < 0)
    {
        nEdges_ = edges().size();
    }

    return nEdges_;
}


inline Foam::label Foam::primitiveMesh::nInternalFaces() const
{
    return nInternalFaces_;
}


inline Foam::label Foam::primitiveMesh::nBoundaryFaces() const
{
    return (nFaces_ - nInternalFaces_);
}


inline Foam::label Foam::primitiveMesh::nFaces() const
{
    return nFaces_;
}


inline Foam::label Foam::primitiveMesh::nCells() const
{
    return nCells_;
}


inline bool Foam::primitiveMesh::isInternalFace(const label faceIndex) const
{
    return faceIndex < nInternalFaces_;
}


inline bool Foam::primitiveMesh::hasCellShapes() const
{
    return cellShapesPtr_;
}


inline bool Foam::primitiveMesh::hasEdges() const
{
    return edgesPtr_;
}


inline bool Foam::primitiveMesh::hasCellCells() const
{
    return ccPtr_;
}


inline bool Foam::primitiveMesh::hasEdgeCells() const
{
    return ecPtr_;
}


inline bool Foam::primitiveMesh::hasPointCells() const
{
    return pcPtr_;
}


inline bool Foam::primitiveMesh::hasCells() const
{
    return cfPtr_;
}


inline bool Foam::primitiveMesh::hasEdgeFaces() const
{
    return efPtr_;
}


inline bool Foam::primitiveMesh::hasPointFaces() const
{
    return pfPtr_;
}


inline bool Foam::primitiveMesh::hasCellEdges() const
{
    return cePtr_;
}


inline bool Foam::primitiveMesh::hasFaceEdges() const
{
    return fePtr_;
}


inline bool Foam::primitiveMesh::hasPointEdges() const
{
    return pePtr_;
}


inline bool Foam::primitiveMesh::hasPointPoints() const
{
    return ppPtr_;
}


inline bool Foam::primitiveMesh::hasCellPoints() const
{
    return cpPtr_;
}


inline bool Foam::primitiveMesh::hasCellCentres() const
{
    return cellCentresPtr_;
}


inline bool Foam::primitiveMesh::hasFaceCentres() const
{
    return faceCentresPtr_;
}


inline bool Foam::primitiveMesh::hasCellVolumes() const
{
    return cellVolumesPtr_;
}


inline bool Foam::primitiveMesh::hasFaceAreas() const
{
    return faceAreasPtr_;
}


// ************************************************************************* //
