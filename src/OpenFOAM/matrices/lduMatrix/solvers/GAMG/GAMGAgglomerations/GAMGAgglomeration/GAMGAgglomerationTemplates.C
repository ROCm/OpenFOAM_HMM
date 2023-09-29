/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "GAMGAgglomeration.H"
#include "mapDistribute.H"
#include "globalIndex.H"


#ifndef OMP_UNIFIED_MEMORY_REQUIRED
#pragma omp requires unified_shared_memory
#define OMP_UNIFIED_MEMORY_REQUIRED
#endif

#define USM_GAMGA_1  

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::GAMGAgglomeration::restrictField
(
    Field<Type>& cf,
    const Field<Type>& ff,
    const labelList& fineToCoarse
) const
{
    cf = Zero;
    label loop_len = ff.size();
    //forAll(ff, i) // macro: for (label i = 0; i < ff.size(); ++i)
    #pragma omp target teams distribute parallel for if(target:loop_len>10000)
    for (label i = 0; i < loop_len; ++i)
    {
       #pragma omp atomic    
        cf[fineToCoarse[i]] += ff[i];
    }
}


template<class Type>
void Foam::GAMGAgglomeration::restrictField
(
    Field<Type>& cf,
    const Field<Type>& ff,
    const label fineLevelIndex,
    const bool procAgglom
) const
{
    const labelList& fineToCoarse = restrictAddressing_[fineLevelIndex];

    if (!procAgglom && ff.size() != fineToCoarse.size())
    {
        FatalErrorInFunction
            << "field does not correspond to level " << fineLevelIndex
            << " sizes: field = " << ff.size()
            << " level = " << fineToCoarse.size()
            << abort(FatalError);
    }

    restrictField(cf, ff, fineToCoarse);

    const label coarseLevelIndex = fineLevelIndex+1;

    if (procAgglom && hasProcMesh(coarseLevelIndex))
    {
        const label coarseComm =
            UPstream::parent(procCommunicator_[coarseLevelIndex]);

        const List<label>& procIDs = agglomProcIDs(coarseLevelIndex);
        const labelList& offsets = cellOffsets(coarseLevelIndex);

        globalIndex::gather
        (
            offsets,
            coarseComm,
            procIDs,
            cf,
            UPstream::msgType(),
            Pstream::commsTypes::nonBlocking    //Pstream::commsTypes::scheduled
        );
    }
}


template<class Type>
void Foam::GAMGAgglomeration::restrictFaceField
(
    Field<Type>& cf,
    const Field<Type>& ff,
    const label fineLevelIndex
) const
{
    const labelList& fineToCoarse = faceRestrictAddressing_[fineLevelIndex];

    if (ff.size() != fineToCoarse.size())
    {
        FatalErrorInFunction
            << "field does not correspond to level " << fineLevelIndex
            << " sizes: field = " << ff.size()
            << " level = " << fineToCoarse.size()
            << abort(FatalError);
    }

    cf = Zero;

    label loop_len = fineToCoarse.size();

    //forAll(fineToCoarse, ffacei)
    #pragma omp target teams distribute parallel for if(target:loop_len > 10000)
    for (label ffacei=0; ffacei < loop_len; ++ffacei)
    {
        label cFace = fineToCoarse[ffacei];

        if (cFace >= 0)
        {
	    #pragma omp atomic	
            cf[cFace] += ff[ffacei];
        }
    }
}


template<class Type>
void Foam::GAMGAgglomeration::prolongField
(
    Field<Type>& ff,
    const Field<Type>& cf,
    const label levelIndex,
    const bool procAgglom
) const
{
    const labelList& fineToCoarse = restrictAddressing_[levelIndex];

    const label coarseLevelIndex = levelIndex+1;

    if (procAgglom && hasProcMesh(coarseLevelIndex))
    {
        const label coarseComm =
            UPstream::parent(procCommunicator_[coarseLevelIndex]);

        const List<label>& procIDs = agglomProcIDs(coarseLevelIndex);
        const labelList& offsets = cellOffsets(coarseLevelIndex);

        const label localSize = nCells_[levelIndex];

        Field<Type> allCf(localSize);
        globalIndex::scatter
        (
            offsets,
            coarseComm,
            procIDs,
            cf,
            allCf,
            UPstream::msgType(),
            Pstream::commsTypes::nonBlocking    //Pstream::commsTypes::scheduled
        );

        

        #pragma omp target teams distribute parallel for if(target:fineToCoarse.size() > 10000)
        for (label i = 0; i < fineToCoarse.size(); ++i)
        //forAll(fineToCoarse, i)
        {
            ff[i] = allCf[fineToCoarse[i]];
        }
    }
    else
    {
	#pragma omp target teams distribute parallel for if(target:fineToCoarse.size() > 10000)
        for (label i = 0; i < fineToCoarse.size(); ++i)
        //forAll(fineToCoarse, i)
        {
            ff[i] = cf[fineToCoarse[i]];
        }
    }
}


// ************************************************************************* //
