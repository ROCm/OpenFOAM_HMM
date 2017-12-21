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

#include "foamVtkAdaptors.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(foamPvCore, 0);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::Ostream& Foam::foamPvCore::printDataArraySelection
(
    Ostream& os,
    vtkDataArraySelection* select
)
{
    if (!select)
    {
        return os;
    }

    const int n = select->GetNumberOfArrays();

    os << n << '(';
    for (int i=0; i < n; ++i)
    {
        if (i) os << ' ';
        os  << select->GetArrayName(i) << '='
            << (select->GetArraySetting(i) ? 1 : 0);
    }
    os << ')';

    return os;
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

    // Name the output block when assigning dataset 0
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
Foam::foamPvCore::getSelectedArraySet
(
    vtkDataArraySelection* select
)
{
    const int n = select->GetNumberOfArrays();
    HashSet<string> enabled(2*n);

    for (int i=0; i < n; ++i)
    {
        if (select->GetArraySetting(i))
        {
            enabled.insert(select->GetArrayName(i));
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

        for (auto k : enabled)
        {
            Info<< " " << k;
        }
        Info<< " )\n";
    }

    return enabled;
}


Foam::Map<Foam::string>
Foam::foamPvCore::getSelectedArrayMap
(
    vtkDataArraySelection* select
)
{
    const int n = select->GetNumberOfArrays();
    Map<string> enabled(2*n);

    for (int i=0; i < n; ++i)
    {
        if (select->GetArraySetting(i))
        {
            enabled.insert(i, select->GetArrayName(i));
        }
    }

    return enabled;
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

    return word::null;
}


void Foam::foamPvCore::printMemory()
{
    memInfo mem;

    if (mem.valid())
    {
        Info<< "mem peak/size/rss: " << mem << endl;
    }
}


vtkSmartPointer<vtkCellArray> Foam::foamPvCore::identityVertices
(
    const label size
)
{
    // VTK_VERTEX
    auto cells = vtkSmartPointer<vtkCellArray>::New();

    UList<vtkIdType> cellsUL = vtkUList(cells, size, 2*size);

    // Cell connectivity for vertex
    // [size, ids.., size, ids...]
    // which means
    // [1, id, 1, id, ...]
    label idx = 0;
    for (label id=0; id < size; ++id)
    {
        cellsUL[idx++] = 1;
        cellsUL[idx++] = id;
    }

    return cells;
}


// ************************************************************************* //
