/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "symmetryPlanePolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(symmetryPlanePolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, symmetryPlanePolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, symmetryPlanePolyPatch, dictionary);
}

// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::symmetryPlanePolyPatch::symmetryPlanePolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, size, start, index, bm, patchType),
    n_(vector::zero)
{}


Foam::symmetryPlanePolyPatch::symmetryPlanePolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, dict, index, bm, patchType),
    n_(vector::zero)
{}


Foam::symmetryPlanePolyPatch::symmetryPlanePolyPatch
(
    const symmetryPlanePolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    n_(vector::zero)
{}


Foam::symmetryPlanePolyPatch::symmetryPlanePolyPatch
(
    const symmetryPlanePolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart),
    n_(vector::zero)
{}


Foam::symmetryPlanePolyPatch::symmetryPlanePolyPatch
(
    const symmetryPlanePolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    polyPatch(pp, bm, index, mapAddressing, newStart),
    n_(vector::zero)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::vector& Foam::symmetryPlanePolyPatch::n() const
{
    // If the symmetry normal is not set calculate it
    // as the average face-normal
    if (magSqr(n_) < 0.5)
    {
        const vectorField& nf(faceNormals());
        n_ = gAverage(nf);

        // Check the symmetry plane is planar
        forAll(nf, facei)
        {
            if (magSqr(n_ - nf[facei]) > SMALL)
            {
                FatalErrorIn("symmetryPlanePolyPatch::n()")
                    << "Symmetry plane '" << name() << "' is not planar"
                    << exit(FatalError);
            }
        }
    }

    return n_;
}


// ************************************************************************* //
