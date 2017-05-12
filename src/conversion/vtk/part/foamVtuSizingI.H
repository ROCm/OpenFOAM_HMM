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

inline bool Foam::foamVtuSizing::decompose() const
{
    return decompose_;
}


inline Foam::label Foam::foamVtuSizing::nCells() const
{
    return nCells_;
}


inline Foam::label Foam::foamVtuSizing::nPoints() const
{
    return nPoints_;
}


inline Foam::label Foam::foamVtuSizing::nVertLabels() const
{
    return nVertLabels_;
}


inline Foam::label Foam::foamVtuSizing::nFaceLabels() const
{
    return nFaceLabels_;
}


inline Foam::label Foam::foamVtuSizing::nCellsPoly() const
{
    return nCellsPoly_;
}


inline Foam::label Foam::foamVtuSizing::nVertPoly() const
{
    return nVertPoly_;
}


inline Foam::label Foam::foamVtuSizing::nAddCells() const
{
    return nAddCells_;
}


inline Foam::label Foam::foamVtuSizing::nAddPoints() const
{
    return nAddPoints_;
}


inline Foam::label Foam::foamVtuSizing::nAddVerts() const
{
    return nAddVerts_;
}


inline Foam::label Foam::foamVtuSizing::nFieldCells() const
{
    return nCells_ + nAddCells_;
}


inline Foam::label Foam::foamVtuSizing::nFieldPoints() const
{
    return nPoints_ + nAddPoints_;
}


inline Foam::label Foam::foamVtuSizing::sizeLegacy() const
{
    return slotSize(contentType::LEGACY, slotType::CELLS);
}


inline Foam::label Foam::foamVtuSizing::sizeLegacy
(
    const enum slotType slot
) const
{
    return slotSize(contentType::LEGACY, slot);
}


inline Foam::label Foam::foamVtuSizing::sizeXml
(
    const enum slotType slot
) const
{
    return slotSize(contentType::XML, slot);
}


inline Foam::label Foam::foamVtuSizing::sizeInternal
(
    const enum slotType slot
) const
{
    return slotSize(contentType::INTERNAL, slot);
}


// ************************************************************************* //
