/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "foamPvCore.H"

#include "memInfo.H"
#include "DynamicList.H"

#include "vtkDataArraySelection.h"
#include "vtkDataSet.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkInformation.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(foamPvCore, 0);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::foamPvCore::addToBlock
(
    vtkMultiBlockDataSet* output,
    vtkDataSet* dataset,
    const arrayRange& selector,
    const label datasetNo,
    const std::string& datasetName
)
{
    const int blockNo = selector.block();

    vtkDataObject* dataObj = output->GetBlock(blockNo);
    vtkMultiBlockDataSet* block = vtkMultiBlockDataSet::SafeDownCast(dataObj);

    if (!block)
    {
        if (dataObj)
        {
            FatalErrorInFunction
                << "Block already has a vtkDataSet assigned to it"
                << endl;
            return;
        }

        block = vtkMultiBlockDataSet::New();
        output->SetBlock(blockNo, block);
        block->Delete();
    }

    if (debug)
    {
        Info<< "block[" << blockNo << "] has "
            << block->GetNumberOfBlocks()
            <<  " datasets prior to adding set " << datasetNo
            <<  " with name: " << datasetName << endl;
    }

    block->SetBlock(datasetNo, dataset);

    // name the output block when assigning dataset 0
    if (datasetNo == 0)
    {
        output->GetMetaData(blockNo)->Set
        (
            vtkCompositeDataSet::NAME(),
            selector.name()
        );
    }

    if (datasetName.size())
    {
        block->GetMetaData(datasetNo)->Set
        (
            vtkCompositeDataSet::NAME(),
            datasetName.c_str()
        );
    }
}


int Foam::foamPvCore::getNumberOfDataSets
(
    vtkMultiBlockDataSet* output,
    const arrayRange& selector
)
{
    const int blockNo = selector.block();

    vtkMultiBlockDataSet* block = vtkMultiBlockDataSet::SafeDownCast
    (
        output->GetBlock(blockNo)
    );

    if (block)
    {
        return block->GetNumberOfBlocks();
    }

    return 0;
}


int Foam::foamPvCore::getSelected
(
    boolList& status,
    vtkDataArraySelection* selection
)
{
    const int n = selection->GetNumberOfArrays();
    if (status.size() != n)
    {
        status.setSize(n);
        status = false;
    }

    int count = 0;
    forAll(status, i)
    {
        const bool setting = selection->GetArraySetting(i);
        if (setting)
        {
            ++count;
        }
        status[i] = setting;
    }

    return count;
}


Foam::hashedWordList Foam::foamPvCore::getSelected
(
    vtkDataArraySelection* select
)
{
    const int n = select->GetNumberOfArrays();
    DynamicList<word> selected(n);

    for (int i=0; i < n; ++i)
    {
        if (select->GetArraySetting(i))
        {
            selected.append(getFirstWord(select->GetArrayName(i)));
        }
    }

    return hashedWordList(selected, true);
}


Foam::hashedWordList Foam::foamPvCore::getSelected
(
    vtkDataArraySelection* select,
    const arrayRange& selector
)
{
    const int n = select->GetNumberOfArrays();
    DynamicList<word> selected(n);

    for (int i = selector.start(); i < selector.end(); ++i)
    {
        if (select->GetArraySetting(i))
        {
            selected.append(getFirstWord(select->GetArrayName(i)));
        }
    }

    return hashedWordList(selected, true);
}


Foam::stringList Foam::foamPvCore::getSelectedArrayEntries
(
    vtkDataArraySelection* select
)
{
    stringList selections(select->GetNumberOfArrays());
    label nElem = 0;

    forAll(selections, elemI)
    {
        if (select->GetArraySetting(elemI))
        {
            selections[nElem++] = select->GetArrayName(elemI);
        }
    }
    selections.setSize(nElem);

    if (debug > 1)
    {
        const int n = select->GetNumberOfArrays();
        Info<< "available(";
        for (int i=0; i < n; ++i)
        {
            Info<< " \"" << select->GetArrayName(i) << "\"";
        }
        Info<< " )\nselected(";

        forAll(selections, elemI)
        {
            Info<< " " << selections[elemI];
        }
        Info<< " )\n";
    }

    return selections;
}


Foam::stringList Foam::foamPvCore::getSelectedArrayEntries
(
    vtkDataArraySelection* select,
    const arrayRange& selector
)
{
    stringList selections(selector.size());
    label nElem = 0;

    for (int i = selector.start(); i < selector.end(); ++i)
    {
        if (select->GetArraySetting(i))
        {
            selections[nElem++] = select->GetArrayName(i);
        }
    }
    selections.setSize(nElem);

    if (debug > 1)
    {
        Info<< "available(";
        for (int i = selector.start(); i < selector.end(); ++i)
        {
            Info<< " \"" << select->GetArrayName(i) << "\"";
        }
        Info<< " )\nselected(";

        forAll(selections, elemI)
        {
            Info<< " " << selections[elemI];
        }
        Info<< " )\n";
    }

    return selections;
}


void Foam::foamPvCore::setSelectedArrayEntries
(
    vtkDataArraySelection* select,
    const stringList& selections
)
{
    const int n = select->GetNumberOfArrays();
    select->DisableAllArrays();

    // Loop through entries, setting values from selectedEntries
    for (int i=0; i < n; ++i)
    {
        const string arrayName(select->GetArrayName(i));

        forAll(selections, elemI)
        {
            if (selections[elemI] == arrayName)
            {
                select->EnableArray(arrayName.c_str());
                break;
            }
        }
    }
}


Foam::word Foam::foamPvCore::getFirstWord(const char* str)
{
    if (str)
    {
        label n = 0;
        while (str[n] && word::valid(str[n]))
        {
            ++n;
        }
        // don't need to re-check for invalid chars
        return word(str, n, false);
    }
    else
    {
        return word::null;
    }
}


void Foam::foamPvCore::printMemory()
{
    memInfo mem;

    if (mem.valid())
    {
        Info<< "mem peak/size/rss: " << mem << endl;
    }
}

// ************************************************************************* //
