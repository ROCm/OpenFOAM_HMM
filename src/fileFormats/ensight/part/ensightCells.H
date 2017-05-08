/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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

Class
    Foam::ensightCells

Description
    Sorting/classification of cells (3D) into corresponding ensight element
    types.

\*---------------------------------------------------------------------------*/

#ifndef ensightCells_H
#define ensightCells_H

#include "labelList.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                        Class ensightCells Declaration
\*---------------------------------------------------------------------------*/

class ensightCells
{
public:

    // Public data

        //- Addressable ensight element types
        enum elemType
        {
            TETRA4,     //!< "tetra4"
            PYRAMID5,   //!< "pyramid5"
            PENTA6,     //!< "penta6"
            HEXA8,      //!< "hexa8"
            NFACED      //!< "nfaced"
        };

        //- Number of element types (5)
        static const label nTypes;

        //- The ensight element type names
        static const char* elemNames[5];


    // Static Member Functions

        //- Return the ensight element name for the specified type
        inline static const char* key(const enum elemType);


private:

    // Private data

        //- Location within a list.
        //  The ensight part number is typically this value +1.
        label index_;

        //- Linear list of ids, sub-sectioned per element type by sub-lists
        labelList address_;

        //- Slices (sub-lists) of the address ids for each element type.
        FixedList<labelRange, 5> slices_;

        //- List of global sizes for each element type.
        //  Used temporarily for local sizes when building the element lists.
        FixedList<label, 5> sizes_;


    // Private Member Functions

        //- Use temporarily stored sizes to redimension the element lists
        void resizeAll();

        //- Disallow default bitwise assignment
        void operator=(const ensightCells&) = delete;


public:

    // Constructors

        //- Construct null, optionally with index
        ensightCells(label partIndex = 0);

        //- Copy constructor. Needed for lists etc.
        ensightCells(const ensightCells& obj);


    //- Destructor
    ~ensightCells();


    // Member Functions

      // Access

        //- The index in a list.
        inline label index() const;

        //- The index in a list, non-const access.
        inline label& index();

        //- The processor local size of all elements.
        inline label size() const;

        //- The processor local size of the specified element type.
        inline label size(const enum elemType) const;

        //- The global number of the specified element type.
        //  This value is only meaningful after a reduce operation.
        inline label total(const enum elemType) const;

        //- The global number of all element types.
        //  This value is only meaningful after a reduce operation.
        label total() const;

        //- The global numbers per element type.
        //  This value is only meaningful after a reduce operation.
        inline const FixedList<label, 5>& totals() const;

        //- The processor local sizes per element type.
        FixedList<label, 5> sizes() const;

        //- Processor local starting offset of element type.
        inline label offset(const enum elemType what) const;

        //- Return the (local) cell ids of the specified element type
        inline const labelUList cellIds(const enum elemType) const;

        //- Return the cell ids of all elements
        inline const labelUList& cellIds() const;


      // Edit

        //- Classify cell types and set the element lists.
        //  The optional indirect addressing can be used when classifying
        //  groups of cells (eg, from a cellZone etc).
        void classify
        (
            const polyMesh& mesh,
            const labelUList& addressing = labelUList::null()
        );


        //- Set addressable sizes to zero, free up addressing memory.
        void clear();

        //- Sum element counts across all processes.
        void reduce();

        //- Sort element lists numerically.
        void sort();


    // Member operators

        //- Return id from linear list of addressing.
        inline label operator[](const label i) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ensightCellsI.H"

#endif

// ************************************************************************* //
