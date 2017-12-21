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

#include "semiImplicitOversetGAMGInterfaceField.H"
#include "addToRunTimeSelectionTable.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(semiImplicitOversetGAMGInterfaceField, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        semiImplicitOversetGAMGInterfaceField,
        lduInterfaceField
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::semiImplicitOversetGAMGInterfaceField::
semiImplicitOversetGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const lduInterfaceField& fineInterface
)
:
    oversetGAMGInterfaceField(GAMGCp, fineInterface)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::semiImplicitOversetGAMGInterfaceField::
~semiImplicitOversetGAMGInterfaceField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::semiImplicitOversetGAMGInterfaceField::updateInterfaceMatrix
(
    scalarField& result,
    const bool add,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    // Add remote values
    if (oversetInterface_.master())
    {
        const labelListList& stencil = oversetInterface_.stencil();

        if (stencil.size() != psiInternal.size())
        {
            FatalErrorInFunction << "psiInternal:" << psiInternal.size()
                << " stencil:" << stencil.size() << exit(FatalError);
        }

        const mapDistribute& map = oversetInterface_.cellInterpolationMap();
        const List<scalarList>& wghts =
            oversetInterface_.cellInterpolationWeights();
        const labelList& cellIDs = oversetInterface_.interpolationCells();
        const scalarList& factor = oversetInterface_.cellInterpolationWeight();
        const scalarField& normalisation = oversetInterface_.normalisation();

        scalarField work(psiInternal);
        map.mapDistributeBase::distribute(work, UPstream::msgType()+1);

        forAll(cellIDs, i)
        {
            label celli = cellIDs[i];

            const scalarList& w = wghts[celli];
            const labelList& nbrs = stencil[celli];
            const scalar f = factor[celli];

            scalar s(0.0);
            forAll(nbrs, nbrI)
            {
                label slotI = nbrs[nbrI];
                s += w[nbrI]*work[slotI];
            }
            if (add)
            {
                s = -1.0*s;
            }

            result[celli] += normalisation[celli]*f*s;
        }
    }
}


// ************************************************************************* //
