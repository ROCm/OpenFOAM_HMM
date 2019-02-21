/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::ijkMesh::ijkMesh()
:
    ijkAddressing()
{}


inline Foam::ijkMesh::ijkMesh(const labelVector& ijk)
:
    ijkAddressing(ijk)
{}


inline Foam::ijkMesh::ijkMesh
(
    const label nx,
    const label ny,
    const label nz
)
:
    ijkAddressing(nx, ny, nz)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::ijkMesh::nPoints() const
{
    if (ijkAddressing::empty())
    {
        return 0;
    }

    const labelVector& n = ijkAddressing::sizes();

    return ((n.x()+1) * (n.y()+1) * (n.z()+1));
}


inline Foam::label Foam::ijkMesh::nCells() const
{
    return ijkAddressing::size();
}


inline Foam::label Foam::ijkMesh::nFaces() const
{
    if (ijkAddressing::empty())
    {
        return 0;
    }

    const labelVector& n = ijkAddressing::sizes();

    return
    (
        ((n.x()+1) * n.y() * n.z())
      + ((n.y()+1) * n.z() * n.x())
      + ((n.z()+1) * n.x() * n.y())
    );
}


inline Foam::label Foam::ijkMesh::nInternalFaces() const
{
    if (ijkAddressing::empty())
    {
        return 0;
    }

    const labelVector& n = ijkAddressing::sizes();

    return
    (
        ((n.x()-1) * n.y() * n.z())
      + ((n.y()-1) * n.z() * n.x())
      + ((n.z()-1) * n.x() * n.y())
    );
}


inline Foam::label Foam::ijkMesh::nBoundaryFaces() const
{
    if (ijkAddressing::empty())
    {
        return 0;
    }

    const labelVector& n = ijkAddressing::sizes();

    return
    (
        (2 * n.y() * n.z())
      + (2 * n.z() * n.x())
      + (2 * n.x() * n.y())
    );
}


inline Foam::label Foam::ijkMesh::cellLabel
(
    const label i,
    const label j,
    const label k
) const
{
    return ijkAddressing::index(i, j, k);
}


inline Foam::label Foam::ijkMesh::pointLabel
(
    const label i,
    const label j,
    const label k
) const
{
    #ifdef FULLDEBUG
    checkIndex(i, j, k, true);
    #endif

    const labelVector& n = sizes();

    return (i + ((n.x()+1) * (j + (n.y()+1) * k)));
}


// ************************************************************************* //
