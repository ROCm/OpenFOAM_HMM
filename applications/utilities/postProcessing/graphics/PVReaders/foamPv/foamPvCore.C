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
#include "vtkSmartPointer.h"

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
    vtkSmartPointer<vtkMultiBlockDataSet> block =
        vtkMultiBlockDataSet::SafeDownCast(dataObj);

    if (!block)
    {
        if (dataObj)
        {
            FatalErrorInFunction
                << "Block already has a vtkDataSet assigned to it"
                << endl;
            return;
        }

        block = vtkSmartPointer<vtkMultiBlockDataSet>::New();
        output->SetBlock(blockNo, block);
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
            selected.append(getFoamName(select->GetArrayName(i)));
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

    for (auto i : selector)
    {
        if (select->GetArraySetting(i))
        {
            selected.append(getFoamName(select->GetArrayName(i)));
        }
    }

    return hashedWordList(selected, true);
}


Foam::HashSet<Foam::string>
Foam::foamPvCore::getSelectedArrayEntries
(
    vtkDataArraySelection* select
)
{
    const int n = select->GetNumberOfArrays();
    HashSet<string> selections(2*n);

    for (int i=0; i < n; ++i)
    {
        if (select->GetArraySetting(i))
        {
            selections.insert(select->GetArrayName(i));
        }
    }

    if (debug > 1)
    {
        const int n = select->GetNumberOfArrays();
        Info<< "available(";
        for (int i=0; i < n; ++i)
        {
            Info<< " \"" << select->GetArrayName(i) << "\"";
        }
        Info<< " )\nselected(";

        for (auto k : selections)
        {
            Info<< " " << k;
        }
        Info<< " )\n";
    }

    return selections;
}


Foam::HashSet<Foam::string>
Foam::foamPvCore::getSelectedArrayEntries
(
    vtkDataArraySelection* select,
    const arrayRange& slice
)
{
    const int n = select->GetNumberOfArrays();
    HashSet<string> enabled(2*n);

    for (auto i : slice)
    {
        if (select->GetArraySetting(i))
        {
            enabled.insert(select->GetArrayName(i));
        }
    }

    if (debug > 1)
    {
        Info<< "available(";
        for (auto i : slice)
        {
            Info<< " \"" << select->GetArrayName(i) << "\"";
        }
        Info<< " )\nselected(";

        for (auto k : enabled)
        {
            Info<< " " << k;
        }
        Info<< " )\n";
    }

    return enabled;
}


void Foam::foamPvCore::setSelectedArrayEntries
(
    vtkDataArraySelection* select,
    const HashSet<string>& enabled
)
{
    const int n = select->GetNumberOfArrays();
    // disable everything not explicitly enabled
    select->DisableAllArrays();

    // Loop through entries, enabling as required
    for (int i=0; i < n; ++i)
    {
        const char* arrayName = select->GetArrayName(i);
        if (enabled.found(arrayName))
        {
            select->EnableArray(arrayName);
        }
    }
}


Foam::word Foam::foamPvCore::getFoamName(const std::string& str)
{
    if (str.size())
    {
        std::string::size_type beg = str.rfind('/');
        if (beg == std::string::npos)
        {
            beg = 0;
        }
        else
        {
            ++beg;
        }

        std::string::size_type end = beg;

        while (str[end] && word::valid(str[end]))
        {
            ++end;
        }

        // Already checked for valid/invalid chars
        return word(str.substr(beg, beg+end), false);
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
