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
#include "symmetryPolyPatch.H"

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
    n_(calcNormal())
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
    n_(calcNormal())
{}


Foam::symmetryPlanePolyPatch::symmetryPlanePolyPatch
(
    const symmetryPlanePolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    n_(calcNormal())
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
    n_(calcNormal())
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
    n_(calcNormal())
{}


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::symmetryPlanePolyPatch::calcNormal() const
{
    if (returnReduce(size(), sumOp<label>()) == 0)
    {
        // No faces in patch. Avoid gAverage complaining and set
        // normal to nonsense value to catch any use
        return vector::rootMax;
    }
    else
    {
        const vectorField& nf(faceNormals());
        vector n = gAverage(nf);

        // Check the symmetry plane is planar
        forAll(nf, facei)
        {
            if (magSqr(n - nf[facei]) > SMALL)
            {
                FatalErrorIn("symmetryPlanePolyPatch::n()")
                    << "Symmetry plane '" << name() << "' is not planar."
                    << endl
                    << " Either split the patch into planar parts"
                    << " or use the " << symmetryPolyPatch::typeName
                    << " patch type"
                    << exit(FatalError);
            }
        }

        return n;
    }
}


// ************************************************************************* //
