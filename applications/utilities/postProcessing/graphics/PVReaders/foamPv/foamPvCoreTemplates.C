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

#include "IOobjectList.H"
#include "vtkDataArraySelection.h"
#include "vtkMultiBlockDataSet.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Type* Foam::foamPvCore::getDataFromBlock
(
    vtkMultiBlockDataSet* output,
    const arrayRange& selector,
    const label datasetNo
)
{
    const int blockNo = selector.block();

    vtkMultiBlockDataSet* block =
    (
        datasetNo < 0
      ? nullptr
      : vtkMultiBlockDataSet::SafeDownCast(output->GetBlock(blockNo))
    );

    if (block)
    {
        return Type::SafeDownCast(block->GetBlock(datasetNo));
    }

    return nullptr;
}


template<class Type>
Foam::label Foam::foamPvCore::addToSelection
(
    vtkDataArraySelection *select,
    const IOobjectList& objects,
    const string& suffix
)
{
    const wordList names = objects.sortedNames(Type::typeName);

    forAll(names, i)
    {
        if (suffix.empty())
        {
            select->AddArray(names[i].c_str());
        }
        else
        {
            select->AddArray((names[i] + suffix).c_str());
        }
    }

    return names.size();
}


// ************************************************************************* //
