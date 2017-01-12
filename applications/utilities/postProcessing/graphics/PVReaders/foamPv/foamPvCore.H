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

Description
    Helpers for OpenFOAM reader interfaces in ParaView.

SourceFiles
    foamPvCore.C

\*---------------------------------------------------------------------------*/

#ifndef foamPvCore_H
#define foamPvCore_H

#include "className.H"
#include "fileName.H"
#include "stringList.H"
#include "boolList.H"
#include "pointList.H"
#include "wordList.H"
#include "Hash.H"
#include "hashedWordList.H"

#include "vtkPoints.h"

// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

class vtkDataArraySelection;
class vtkDataSet;
class vtkMultiBlockDataSet;
class vtkIndent;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class IOobjectList;

/*---------------------------------------------------------------------------*\
                         Class foamPvCore Declaration
\*---------------------------------------------------------------------------*/

class foamPvCore
{
public:

    //- Bookkeeping for GUI checklists and multi-block organization
    //  Works like a SubList selection.
    class arrayRange
    {
        const char *name_;
        int block_;
        int start_;
        int size_;

    public:

        //- Construct with given name for the specified block
        arrayRange(const char *name, const int blockNo=0)
        :
            name_(name),
            block_(blockNo),
            start_(0),
            size_(0)
        {}

        //- Return the block holding these datasets
        int block() const
        {
            return block_;
        }

        //- Assign block number, return previous value
        int block(int blockNo)
        {
            int prev = block_;
            block_ = blockNo;
            return prev;
        }

        //- Return the name
        const char* name() const
        {
            return name_;
        }

        //- The array start index
        int start() const
        {
            return start_;
        }

        //- The array end index
        int end() const
        {
            return start_ + size_;
        }

        //- The sub-list size
        int size() const
        {
            return size_;
        }

        //- True if the sub-list is empty
        bool empty() const
        {
            return !size_;
        }

        //- Reset the size to zero and optionally assign a new start
        void reset(const int startAt = 0)
        {
            start_ = startAt;
            size_  = 0;
        }

        //- Assign new start and reset the size
        void operator=(const int i)
        {
            reset(i);
        }

        //- Increment the size
        void operator+=(const int n)
        {
            size_ += n;
        }

    }; // End class arrayRange


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        foamPvCore(const foamPvCore&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const foamPvCore&) = delete;

public:

    //- Static data members
    ClassName("foamPvCore");


    //- Construct null
    foamPvCore()
    {}


    //- Convenience method for the VTK multiblock API
    static void addToBlock
    (
        vtkMultiBlockDataSet* output,
        vtkDataSet* dataset,
        const arrayRange& selector,
        const label datasetNo,
        const std::string& datasetName
    );

    //- Convenience method for the VTK multiblock API
    //  Always returns a nullptr if datasetNo is negative
    template<class Type=vtkDataSet>
    static Type* getDataFromBlock
    (
        vtkMultiBlockDataSet* output,
        const arrayRange& selector,
        const label datasetNo
    );

    //- Convenience method for the VTK multiblock API
    static int getNumberOfDataSets
    (
        vtkMultiBlockDataSet* output,
        const arrayRange& selector
    );

    //- Add objects of Type to array selection
    template<class Type>
    static label addToSelection
    (
        vtkDataArraySelection* select,
        const IOobjectList& objects,
        const string& suffix = string::null
    );


    //- Retrieve the current selections into a boolList
    static int getSelected
    (
        boolList& lst,
        vtkDataArraySelection* select
    );

    //- Retrieve the current selections as a wordHashSet
    static hashedWordList getSelected
    (
        vtkDataArraySelection* select
    );


    //- Retrieve a sub-list of the current selections
    static hashedWordList getSelected
    (
        vtkDataArraySelection* select,
        const arrayRange& selector
    );


    //- Retrieve the current selections
    static stringList getSelectedArrayEntries
    (
        vtkDataArraySelection* select
    );

    //- Retrieve a sub-list of the current selections
    static stringList getSelectedArrayEntries
    (
        vtkDataArraySelection* select,
        const arrayRange& selector
    );


    //- Set selection(s)
    static void setSelectedArrayEntries
    (
        vtkDataArraySelection* select,
        const stringList& selections
    );


    //- Extract up to the first non-word characters
    static word getFirstWord(const char* str);

    //- Simple memory used debugging information
    static void printMemory();


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}; // End class foamPvCore

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "foamPvCoreTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
