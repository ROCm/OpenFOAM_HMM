/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

inline const Foam::fvMesh& Foam::fvMeshSubset::baseMesh() const noexcept
{
    return baseMesh_;
}


inline const Foam::fvMesh& Foam::fvMeshSubset::mesh() const noexcept
{
    return fvMeshSubsetPtr_ ? *fvMeshSubsetPtr_ : baseMesh_;
}


inline bool Foam::fvMeshSubset::hasSubMesh() const noexcept
{
    return bool(fvMeshSubsetPtr_);
}


inline const Foam::fvMesh& Foam::fvMeshSubset::subMesh() const
{
    checkCellSubset();

    return *fvMeshSubsetPtr_;
}


inline Foam::fvMesh& Foam::fvMeshSubset::subMesh()
{
    checkCellSubset();

    return *fvMeshSubsetPtr_;
}


inline const Foam::labelList& Foam::fvMeshSubset::pointMap() const
{
    checkCellSubset();

    return pointMap_;
}


inline const Foam::labelList& Foam::fvMeshSubset::faceMap() const
{
    checkCellSubset();

    return faceMap_;
}


inline const Foam::labelList& Foam::fvMeshSubset::faceFlipMap() const
{
    if (!faceFlipMapPtr_)
    {
        calcFaceFlipMap();
    }

    return *faceFlipMapPtr_;
}


inline const Foam::labelList& Foam::fvMeshSubset::cellMap() const
{
    checkCellSubset();

    return cellMap_;
}


inline const Foam::labelList& Foam::fvMeshSubset::patchMap() const
{
    checkCellSubset();

    return patchMap_;
}


// ************************************************************************* //
