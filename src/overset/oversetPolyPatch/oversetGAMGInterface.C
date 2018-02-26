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

#include "oversetGAMGInterface.H"
#include "addToRunTimeSelectionTable.H"
#include "Map.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(oversetGAMGInterface, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterface,
        oversetGAMGInterface,
        lduInterface
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oversetGAMGInterface::oversetGAMGInterface
(
    const label index,
    const lduInterfacePtrsList& coarseInterfaces,
    const lduInterface& fineInterface,
    const labelField& localRestrictAddressing,
    const labelField& neighbourRestrictAddressing,
    const label fineLevelIndex,
    const label coarseComm
)
:
    GAMGInterface(index, coarseInterfaces),
    fineOversetInterface_
    (
        refCast<const oversetLduInterface>(fineInterface)
    )
{
    //Pout<< "Constructing oversetGAMGInterface index:" << index
    //    << " from:" << fineOversetInterface_.name()
    //    << " size:" << localRestrictAddressing.size()
    //    << endl;

    // From coarse face to coarse cell
    DynamicList<label> dynFaceCells(localRestrictAddressing.size());
    // From fine face to coarse face
    DynamicList<label> dynFaceRestrictAddressing(dynFaceCells.size());


    // From coarse cell to coarse face
    Map<label> cellToCoarseFace(dynFaceCells.size());

    // Construct face agglomeration from cell agglomeration
    forAll(localRestrictAddressing, ffi)
    {
        // Get coarse cell
        label coarseCelli = localRestrictAddressing[ffi];

        // Do we have coarse face for it?
        Map<label>::iterator iter = cellToCoarseFace.find(coarseCelli);
        if (iter == cellToCoarseFace.end())
        {
            label coarseFacei = dynFaceCells.size();
            cellToCoarseFace.insert(coarseCelli, coarseFacei);
            dynFaceCells.append(coarseCelli);
            dynFaceRestrictAddressing.append(coarseFacei);
        }
        else
        {
            dynFaceRestrictAddressing.append(iter());
        }
    }

    faceCells_.transfer(dynFaceCells);
    faceRestrictAddressing_.transfer(dynFaceRestrictAddressing);

    //Pout<< "** Constructed interface:" << name()
    //    << " of size:" << size() << endl;


    // Determine this level's stencil from the fine stencil
    if (master())
    {
        const mapDistribute& fineMap =
            fineOversetInterface_.cellInterpolationMap();
        const List<scalarList>& fineWghts =
            fineOversetInterface_.cellInterpolationWeights();
        const labelListList& fineStencil = fineOversetInterface_.stencil();
        const labelList& fineCellIDs =
            fineOversetInterface_.interpolationCells();
        const scalarList& fineFactor =
            fineOversetInterface_.cellInterpolationWeight();
        const scalarField& fineNorm = fineOversetInterface_.normalisation();
        const labelList& restrictMap = fineOversetInterface_.restrictMap();


        if (localRestrictAddressing.size() && restrictMap.empty())
        {
            FatalErrorInFunction
                << "The restrict addressing has not been saved on interface "
                << this->name()
                << ". This gets stored when running internalFieldTransfer"
                << exit(FatalError);
        }


        // 'Interpolate' the restrict map so locally we know what the remote
        // cells are going to be.

        label nCoarseCells = max(restrictMap)+1;
        globalIndex globalNumbering(nCoarseCells);

        labelList globalCoarseIDs(restrictMap.size());
        forAll(restrictMap, fineCelli)
        {
            globalCoarseIDs[fineCelli] =
                globalNumbering.toGlobal(restrictMap[fineCelli]);
        }
        fineMap.distribute(globalCoarseIDs);

        //Pout<< this->name()
        //    << " index:" << index
        //    << " size:" << this->size()
        //    << " restrictMap:" << restrictMap.size()
        //    << " fineNumInterpolate:" << fineCellIDs.size()
        //    << " nCoarseCells:" << nCoarseCells
        //    << endl;



        // Accumulate the coarse level stencil

        // Number of fine cells contributing to the coarse cell
        labelList nFineCells(nCoarseCells, 0);

        stencil_.setSize(nCoarseCells);
        cellInterpolationWeights_.setSize(nCoarseCells);
        cellInterpolationWeight_.setSize(nCoarseCells, 0.0);
        normalisation_.setSize(nCoarseCells, 0.0);

        forAll(fineCellIDs, i)
        {
            label fineCelli = fineCellIDs[i];
            label coarseCelli = restrictMap[fineCelli];

            const scalarList& w = fineWghts[fineCelli];
            const labelList& nbrs = fineStencil[fineCelli];

            // Accumulate stencil

            labelList& coarseStencil = stencil_[coarseCelli];
            scalarList& coarseWghts = cellInterpolationWeights_[coarseCelli];

            label sz = coarseStencil.size();
            coarseStencil.setSize(sz+nbrs.size());
            coarseWghts.setSize(coarseStencil.size());

            forAll(nbrs, i)
            {
                coarseStencil[sz] = globalCoarseIDs[nbrs[i]];
                coarseWghts[sz] = w[i];
                sz++;
            }

            // Accumulate weight
            cellInterpolationWeight_[coarseCelli] += fineFactor[fineCelli];
            normalisation_[coarseCelli] += fineNorm[fineCelli];
            nFineCells[coarseCelli]++;
        }

        // Normalise weight (0 .. 1). Take average? Volume average? Min? Max?
        forAll(nFineCells, coarseCelli)
        {
            label nAvg = nFineCells[coarseCelli];
            if (nAvg > 0)
            {
                cellInterpolationWeight_[coarseCelli] /= nAvg;
            }
        }


        // Work the stencil back from global cells into a mapDistribute
        // + slots
        List<Map<label>> compactMap;
        mapPtr_.reset(new mapDistribute(globalNumbering, stencil_, compactMap));


        // Determine the interpolationCells
        label nInterpolate = 0;
        forAll(stencil_, coarseCelli)
        {
            if (stencil_[coarseCelli].size())
            {
                nInterpolate++;
            }
        }

        //Pout<< this->name() << " nInterpolate:" << nInterpolate << endl;

        interpolationCells_.setSize(nInterpolate);
        nInterpolate = 0;
        forAll(stencil_, coarseCelli)
        {
            if (stencil_[coarseCelli].size())
            {
                interpolationCells_[nInterpolate++] = coarseCelli;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::oversetGAMGInterface::~oversetGAMGInterface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::labelField> Foam::oversetGAMGInterface::internalFieldTransfer
(
    const Pstream::commsTypes,
    const labelUList& restrictMap
) const
{
    // Store the restrictMap. This routine gets used for
    // - GAMGAgglomeration      : this is the use we want to intercept.
    // - GAMGProcAgglomeration  : to find out the cell number on the other side
    // - MGridGenGAMGAgglomeration: same
    if (master())
    {
        restrictMap_ = restrictMap;
    }

    return tmp<labelField>::New(restrictMap, faceCells());
}


// ************************************************************************* //
