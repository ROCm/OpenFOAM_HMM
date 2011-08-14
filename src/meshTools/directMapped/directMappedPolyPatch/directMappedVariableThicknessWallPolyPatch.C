/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "directMappedVariableThicknessWallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(directMappedVariableThicknessWallPolyPatch, 0);

    addToRunTimeSelectionTable
    (
        polyPatch,
        directMappedVariableThicknessWallPolyPatch,
        word
    );

    addToRunTimeSelectionTable
    (
        polyPatch,
        directMappedVariableThicknessWallPolyPatch,
        dictionary
    );
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::directMappedVariableThicknessWallPolyPatch::
directMappedVariableThicknessWallPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
:
    directMappedWallPolyPatch(name, size, start, index, bm),
    thickness_(size)
{}


Foam::directMappedVariableThicknessWallPolyPatch::
directMappedVariableThicknessWallPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const word& sampleRegion,
    const directMappedPatchBase::sampleMode mode,
    const word& samplePatch,
    const vectorField& offset,
    const polyBoundaryMesh& bm
)
:
    directMappedWallPolyPatch(name, size, start, index, bm),
    thickness_(size)
{}


Foam::directMappedVariableThicknessWallPolyPatch::
directMappedVariableThicknessWallPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const word& sampleRegion,
    const directMappedPatchBase::sampleMode mode,
    const word& samplePatch,
    const vector& offset,
    const polyBoundaryMesh& bm
)
:
    directMappedWallPolyPatch(name, size, start, index, bm),
    thickness_(size)
{}


Foam::directMappedVariableThicknessWallPolyPatch::
directMappedVariableThicknessWallPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    directMappedWallPolyPatch(name, dict, index, bm),
    thickness_(scalarField("thickness", dict, this->size()))
{}


Foam::directMappedVariableThicknessWallPolyPatch::
directMappedVariableThicknessWallPolyPatch
(
    const directMappedVariableThicknessWallPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    directMappedWallPolyPatch(pp, bm),
    thickness_(pp.thickness_)
{}


Foam::directMappedVariableThicknessWallPolyPatch::
directMappedVariableThicknessWallPolyPatch
(
    const directMappedVariableThicknessWallPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    directMappedWallPolyPatch(pp, bm, index, newSize, newStart),
    thickness_(newSize)
{}


Foam::directMappedVariableThicknessWallPolyPatch::
directMappedVariableThicknessWallPolyPatch
(
    const directMappedVariableThicknessWallPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    directMappedWallPolyPatch(pp, bm, index, mapAddressing, newStart),
    thickness_(pp.size())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::directMappedVariableThicknessWallPolyPatch::
~directMappedVariableThicknessWallPolyPatch()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::directMappedVariableThicknessWallPolyPatch::
write(Foam::Ostream& os) const
{
    os.writeKeyword("thickness") << thickness_ << token::END_STATEMENT << nl;
}

// ************************************************************************* //
