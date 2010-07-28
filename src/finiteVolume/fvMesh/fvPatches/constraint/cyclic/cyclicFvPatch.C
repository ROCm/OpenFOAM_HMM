/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "cyclicFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "transform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(cyclicFvPatch, 0);
addToRunTimeSelectionTable(fvPatch, cyclicFvPatch, polyPatch);


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Make patch weighting factors
void cyclicFvPatch::makeWeights(scalarField& w) const
{
    const cyclicFvPatch& nbrPatch = neighbFvPatch();

    const scalarField& magFa = magSf();
    const scalarField& nbrMagFa = nbrPatch.magSf();

    scalarField deltas = nf() & fvPatch::delta();
    scalarField nbrDeltas = nbrPatch.nf() & nbrPatch.fvPatch::delta();

    forAll(magFa, facei)
    {
        scalar avFa = (magFa[facei] + nbrMagFa[facei])/2.0;

        if (mag(magFa[facei] - nbrMagFa[facei])/avFa > 1e-4)
        {
            FatalErrorIn("cyclicFvPatch::makeWeights(scalarField& w) const")
                << "face " << facei << " areas do not match by "
                << 100*mag(magFa[facei] - nbrMagFa[facei])/avFa
                << "% -- possible face ordering problem"
                << abort(FatalError);
        }

        scalar di = deltas[facei];
        scalar dni = nbrDeltas[facei];

        w[facei] = dni/(di + dni);
    }
}


// Make patch face - neighbour cell distances
void cyclicFvPatch::makeDeltaCoeffs(scalarField& dc) const
{
    //const cyclicPolyPatch& nbrPatch = cyclicPolyPatch_.neighbPatch();
    const cyclicFvPatch& nbrPatch = neighbFvPatch();

    scalarField deltas = nf() & fvPatch::delta();
    scalarField nbrDeltas = nbrPatch.nf() & nbrPatch.fvPatch::delta();

    forAll(deltas, facei)
    {
        scalar di = deltas[facei];
        scalar dni = nbrDeltas[facei];

        dc[facei] = 1.0/(di + dni);
    }
}


// Return delta (P to N) vectors across coupled patch
tmp<vectorField> cyclicFvPatch::delta() const
{
    vectorField patchD = fvPatch::delta();
    vectorField nbrPatchD = neighbFvPatch().fvPatch::delta();

    tmp<vectorField> tpdv(new vectorField(patchD.size()));
    vectorField& pdv = tpdv();

    // To the transformation if necessary
    if (parallel())
    {
        forAll(patchD, facei)
        {
            vector ddi = patchD[facei];
            vector dni = nbrPatchD[facei];

            pdv[facei] = ddi - dni;
        }
    }
    else
    {
        forAll(patchD, facei)
        {
            vector ddi = patchD[facei];
            vector dni = nbrPatchD[facei];

            pdv[facei] = ddi - transform(forwardT()[0], dni);
        }
    }

    return tpdv;
}


tmp<labelField> cyclicFvPatch::interfaceInternalField
(
    const unallocLabelList& internalData
) const
{
    return patchInternalField(internalData);
}


tmp<labelField> cyclicFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& iF
) const
{
    return neighbFvPatch().patchInternalField(iF);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
