/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "processorFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "polyBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(processorFvPatch, 0);
addToRunTimeSelectionTable(fvPatch, processorFvPatch, polyPatch);


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void processorFvPatch::makeWeights(scalarField& w) const
{
    if (Pstream::parRun())
    {
        const processorPolyPatch& pp = procPolyPatch();

//Pout<< name() << " pp.neighbFaceAreas():" << pp.neighbFaceAreas()
//    << endl;
//
//Pout<< name() << " pp.neighbFaceCentres():" << pp.neighbFaceCentres()
//    << endl;
//
//Pout<< name() << " pp.neighbFaceCellCentres():" << pp.neighbFaceCellCentres()
//    << endl;


        // The face normals point in the opposite direction on the other side
        scalarField neighbFaceCentresCn
        (
            (
                pp.neighbFaceAreas()
               /(mag(pp.neighbFaceAreas()) + VSMALL)
            )
          & (
              pp.neighbFaceCentres()
            - pp.neighbFaceCellCentres()
            )
        );

        w = neighbFaceCentresCn/((nf()&fvPatch::delta()) + neighbFaceCentresCn);
    }
    else
    {
        w = 1.0;
    }
//Pout<< name() << " w:" << w
//    << endl;
}


void processorFvPatch::makeDeltaCoeffs(scalarField& dc) const
{
//Pout<< name() << " fvPatch::delta():" << fvPatch::delta()
//    << endl;
//
//Pout<< name() << " nf():" << nf()
//    << endl;
//
//Pout<< name() << " weights():" << weights()
//    << endl;

    if (Pstream::parRun())
    {
        dc = (1.0 - weights())/(nf() & fvPatch::delta());
    }
    else
    {
        dc = 1.0/(nf() & fvPatch::delta());
    }
//Pout<< name() << " dc:" << dc << endl;
}


tmp<vectorField> processorFvPatch::delta() const
{
    tmp<vectorField> deltaFld(fvPatch::delta());

    if (Pstream::parRun())
    {
        // Do the transformation if necessary.

        const processorPolyPatch& pp = procPolyPatch();

        const pointField& nfc = pp.neighbFaceCentres();
        const pointField& ncc = pp.neighbFaceCellCentres();

        forAll(pp.patchIDs(), i)
        {
            const SubField<point> subFc(pp.subSlice(nfc, i));
            const SubField<point> subCc(pp.subSlice(ncc, i));
            SubField<vector> subDelta(pp.subSlice(deltaFld(), i));
            const vectorField& subFld = static_cast<const vectorField&>
            (
                subDelta
            );

            label patchI = pp.patchIDs()[i];

//Pout<< name() << " delta:" << " subFc:" << subFc
//    << " subCc:" << subCc << " subDelta:" << subDelta
//    << endl;

            if (patchI == -1)
            {
                const_cast<vectorField&>(subFld) -= (subFc - subCc);
            }
            else
            {
                const coupledPolyPatch& subPatch =
                    refCast<const coupledPolyPatch>
                    (
                        pp.boundaryMesh()[patchI]
                    );

                if (subPatch.parallel())
                {
                    const_cast<vectorField&>(subFld) -= (subFc - subCc);
                }
                else
                {
                    const_cast<vectorField&>(subFld) -= transform
                    (
                        subPatch.forwardT(),
                        subFc - subCc
                    );
                }
            }
//Pout<< name() << " subDelta:" << subDelta
//    << endl;

        }
    }

    return deltaFld;
}


tmp<labelField> processorFvPatch::interfaceInternalField
(
    const unallocLabelList& internalData
) const
{
    return patchInternalField(internalData);
}


//void processorFvPatch::initTransfer
//(
//    const Pstream::commsTypes commsType,
//    const unallocLabelList& interfaceData
//) const
//{
//    send(commsType, interfaceData);
//}
//
//
//tmp<labelField> processorFvPatch::transfer
//(
//    const Pstream::commsTypes commsType,
//    const unallocLabelList&
//) const
//{
//    return receive<label>(commsType, this->size());
//}


void processorFvPatch::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& iF
) const
{
    send(commsType, patchInternalField(iF)());
}


tmp<labelField> processorFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList&
) const
{
    return receive<label>(commsType, this->size());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
