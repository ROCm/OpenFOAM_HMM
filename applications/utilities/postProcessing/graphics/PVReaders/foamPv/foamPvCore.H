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
#include "HashSet.H"
#include "Map.H"
#include "hashedWordList.H"
#include "labelRange.H"

// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

class vtkCellArray;
class vtkDataArraySelection;
class vtkDataSet;
class vtkIndent;
class vtkMultiBlockDataSet;
class vtkPoints;
template<class T> class vtkSmartPointer;


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
    :
        public labelRange
    {
        const char *name_;
        int block_;

    public:

        //- Construct with given name for the specified block
        arrayRange(const char *name, int blockNo=0)
        :
            labelRange(),
            name_(name),
            block_(blockNo)
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

        //- Reset the start/size directly
        using labelRange::reset;

        //- Reset the size to zero and optionally assign a new start
        void reset(label startAt = 0)
        {
            clear();
            setStart(startAt);
        }

        //- Increment the size
        void operator+=(label n)
        {
            setSize(size() + n);
        }

        //- True if the labelRange intersects any key in the Map
        template<class T>
        bool intersects(const Map<T>& map) const
        {
            for (const label idx : *this)
            {
                if (map.found(idx))
                {
                    return true;
                }
            }

            return false;
        }

        //- The intersection ids with keys in the Map
        template<class T>
        List<label> intersection(const Map<T>& map) const
        {
            List<label> indices(Foam::min(map.size(), this->size()));

            if (indices.size())
            {
                label n = 0;

                for (const label idx : *this)
                {
                    if (map.found(idx))
                    {
                        indices[n++] = idx;
                    }
                }

                indices.setSize(n);
            }

            return indices;
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
    foamPvCore() = default;


    //- Print information about vtkDataArraySelection
    static Ostream& printDataArraySelection
    (
        Ostream& os,
        vtkDataArraySelection* select
    );


    //- Convenience method for the VTK multiblock API
    static void addToBlock
    (
        vtkMultiBlockDataSet* output,
        vtkDataSet* dataset,
        const arrayRange& selector,
        const label datasetNo,
        const std::string& datasetName
    );

    //- Add names to array selection
    template<class StringType>
    static label addToArray
    (
        vtkDataArraySelection* select,
        const std::string& prefix,
        const UList<StringType>& names
    );

    //- Add names to array selection
    template<class StringType>
    static label addToArray
    (
        vtkDataArraySelection* select,
        const UList<StringType>& names,
        const std::string& suffix = string::null
    );

    //- Add objects of Type to array selection
    template<class Type>
    static label addToSelection
    (
        vtkDataArraySelection* select,
        const std::string& prefix,
        const IOobjectList& objects
    );

    //- Add objects of Type to array selection
    template<class Type>
    static label addToSelection
    (
        vtkDataArraySelection* select,
        const IOobjectList& objects,
        const std::string& suffix = string::null
    );

    //- Add objects of Type to array selection
    template<class Type>
    static label addToSelection
    (
        vtkDataArraySelection* select,
        const std::string& prefix,
        const HashTable<wordHashSet>& objects
    );


    //- Add objects of Type to array selection
    template<class Type>
    static label addToSelection
    (
        vtkDataArraySelection* select,
        const HashTable<wordHashSet>& objects,
        const std::string& suffix = string::null
    );


    //- Retrieve the current selections as a hashedWordList,
    //  while stripping off any prefix or suffix
    static hashedWordList getSelected
    (
        vtkDataArraySelection* select
    );


    //- Retrieve a sub-list of the current selections as a hashedWordList,
    //  while stripping off any prefix or suffix
    static hashedWordList getSelected
    (
        vtkDataArraySelection* select,
        const arrayRange& selector
    );


    //- Retrieve the currently enabled selections as hashset
    static HashSet<string> getSelectedArraySet
    (
        vtkDataArraySelection* select
    );

    //- Retrieve the currently enabled selections as id/name map
    static Map<string> getSelectedArrayMap
    (
        vtkDataArraySelection* select
    );

    //- Enable the selection(s)
    template<class AnyValue, class AnyHasher>
    static void setSelectedArrayEntries
    (
        vtkDataArraySelection* select,
        const HashTable<AnyValue, string, AnyHasher>& enabled
    );


    //- Extract the first word characters after a slash
    static word getFoamName(const std::string& str);


    //- Simple memory used debugging information
    static void printMemory();


    //- Return an identity list of VTK_VERTEX
    static vtkSmartPointer<vtkCellArray> identityVertices
    (
        const label size
    );


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
