/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

#include "foamVtuSizing.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::vtk::vtuSizing::decompose() const
{
    return decompose_;
}


inline Foam::label Foam::vtk::vtuSizing::nCells() const
{
    return nCells_;
}


inline Foam::label Foam::vtk::vtuSizing::nPoints() const
{
    return nPoints_;
}


inline Foam::label Foam::vtk::vtuSizing::nVertLabels() const
{
    return nVertLabels_;
}


inline Foam::label Foam::vtk::vtuSizing::nFaceLabels() const
{
    return nFaceLabels_;
}


inline Foam::label Foam::vtk::vtuSizing::nCellsPoly() const
{
    return nCellsPoly_;
}


inline Foam::label Foam::vtk::vtuSizing::nVertPoly() const
{
    return nVertPoly_;
}


inline Foam::label Foam::vtk::vtuSizing::nAddCells() const
{
    return nAddCells_;
}


inline Foam::label Foam::vtk::vtuSizing::nAddPoints() const
{
    return nAddPoints_;
}


inline Foam::label Foam::vtk::vtuSizing::nAddVerts() const
{
    return nAddVerts_;
}


inline Foam::label Foam::vtk::vtuSizing::nFieldCells() const
{
    return nCells_ + nAddCells_;
}


inline Foam::label Foam::vtk::vtuSizing::nFieldPoints() const
{
    return nPoints_ + nAddPoints_;
}


inline Foam::label Foam::vtk::vtuSizing::sizeLegacy() const
{
    return sizeOf(contentType::LEGACY, slotType::CELLS);
}


inline Foam::label Foam::vtk::vtuSizing::sizeLegacy
(
    const enum slotType slot
) const
{
    return sizeOf(contentType::LEGACY, slot);
}


inline Foam::label Foam::vtk::vtuSizing::sizeXml
(
    const enum slotType slot
) const
{
    return sizeOf(contentType::XML, slot);
}


inline Foam::label Foam::vtk::vtuSizing::sizeInternal
(
    const enum slotType slot
) const
{
    return sizeOf(contentType::INTERNAL, slot);
}


// ************************************************************************* //
