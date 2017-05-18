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
    Foam::ensightFaces

Description
    Sorting/classification of faces (2D) into corresponding ensight types

\*---------------------------------------------------------------------------*/

#ifndef ensightFaces_H
#define ensightFaces_H

#include "boolList.H"
#include "labelList.H"
#include "faceList.H"
#include "FixedList.H"
#include "PackedBoolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ensightFaces Declaration
\*---------------------------------------------------------------------------*/

class ensightFaces
{
public:

    // Public data

        //- Addressable ensight element types
        enum elemType
        {
            TRIA3,      //!< "tria3"
            QUAD4,      //!< "quad4"
            NSIDED      //!< "nsided"
        };

        //- Number of element types (3)
        static const label nTypes;

        //- The ensight element type names
        static const char* elemNames[3];


    // Static Member Functions

        //- Return the ensight element name for the specified type
        static inline const char* key(const enum elemType);


private:

    // Private data

        //- Location within a list.
        //  The ensight part number is typically this value +1.
        label index_;

        //- Linear list of ids, sub-sectioned per element type by sub-lists
        labelList address_;

        //- Linear list of face-flips
        boolList flipMap_;

        //- Slices (sub-lists) of the address and flips for each element type.
        FixedList<labelRange, 3> slices_;

        //- List of global sizes for each element type.
        //  Used temporarily for local sizes when building the element lists.
        FixedList<label, 3> sizes_;


    // Private Member Functions

        //- Simple classifier
        inline static elemType whatType(const face& f);

        //- Low-level internal addition routine
        inline void add(const face& f, const label id, const bool flip = false);

        //- Use temporarily stored sizes to redimension the element lists
        void resizeAll();

        //- Disallow default bitwise assignment
        void operator=(const ensightFaces&) = delete;


public:

    // Constructors

        //- Construct null, optionally with index
        ensightFaces(label partIndex = 0);

        //- Copy constructor. Needed for lists etc.
        ensightFaces(const ensightFaces& obj);


    //- Destructor
    ~ensightFaces();


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

        //- The global number of all element types.
        //  This value is only meaningful after a reduce operation.
        label total() const;

        //- The global number of the specified element type.
        //  This value is only meaningful after a reduce operation.
        inline label total(const enum elemType) const;

        //- The global numbers per element type.
        //  This value is only meaningful after a reduce operation.
        inline const FixedList<label, 3>& totals() const;

        //- The processor local sizes per element type.
        FixedList<label, 3> sizes() const;

        //- Processor local starting offset of element type.
        inline label offset(const enum elemType what) const;

        //- Return the (local) face ids of the specified element type
        inline const labelUList faceIds(const enum elemType) const;

        //- Return the processor local face ids of all elements
        inline const labelUList& faceIds() const;

        //- Return the processor local flip-map of all elements
        inline const boolList& flipMap() const;


      // Edit

        //- Classify the face types, set element list.
        void classify(const faceList& faces);


        //- Classify the face types, set element list.
        //  The indirect addressing can be used when classifying groups of
        //  face (eg, from a faceZone etc) with an optional flipMap.
        //  The optional exclude marker can be used to skip faces on particular
        //  boundary types or regions.
        void classify
        (
            const faceList& faces,
            const labelUList& addressing,
            const boolList& flipMap = boolList(),
            const PackedBoolList& exclude = PackedBoolList()
        );


        //- Set addressable sizes to zero, free up addressing memory.
        void clear();

        //- Sum element counts across all processes.
        void reduce();

        //- Sort element lists numerically.
        void sort();


    // Member operators

        //- Return element from linear-list.
        inline label operator[](const label i) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ensightFacesI.H"

#endif

// ************************************************************************* //
