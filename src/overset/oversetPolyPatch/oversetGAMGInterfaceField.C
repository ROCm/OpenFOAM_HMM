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

#include "oversetGAMGInterfaceField.H"
#include "addToRunTimeSelectionTable.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(oversetGAMGInterfaceField, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        oversetGAMGInterfaceField,
        lduInterfaceField
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oversetGAMGInterfaceField::oversetGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const lduInterfaceField& fineInterface
)
:
    GAMGInterfaceField(GAMGCp, fineInterface),
    oversetInterface_(refCast<const oversetGAMGInterface>(GAMGCp))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::oversetGAMGInterfaceField::~oversetGAMGInterfaceField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::oversetGAMGInterfaceField::updateInterfaceMatrix
(
    scalarField& result,
    const bool add,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    //Pout<< "oversetGAMGInterfaceField::updateInterfaceMatrix: at:"
    //    << oversetInterface_.name()
    //    << " nCells:" << result.size() << endl;
    //
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

            bool hasRemote = false;
            forAll(nbrs, nbrI)
            {
                label slotI = nbrs[nbrI];
                if (slotI >= psiInternal.size())
                {
                    hasRemote = true;
                    break;
                }
            }

            if (hasRemote)
            {
                //Pout<< "oversetGAMGInterfaceField : Interpolating " << celli
                //    << " from remote values (if any):" << endl;
                scalar s(0.0);
                forAll(nbrs, nbrI)
                {
                    label slotI = nbrs[nbrI];
                    if (slotI >= psiInternal.size())
                    {
                        //Pout<< "    remote value " << work[slotI]
                        //    << " from slot " << slotI << " with w:" << w[nbrI]
                        //    << endl;
                        s += w[nbrI]*work[slotI];
                    }
                }
                //Pout<< "oversetGAMGInterfaceField : Interpolated value:"
                //    << s << endl;
                //scalar oldPsi = result[celli];
                if (add)
                {
                    s = -1.0*s;
                }

                result[celli] += normalisation[celli]*f*s;
                //Pout<< "oversetGAMGInterfaceField : result was:" << oldPsi
                //    << " now:" << result[celli] << endl;
            }
        }
    }
}


// ************************************************************************* //
