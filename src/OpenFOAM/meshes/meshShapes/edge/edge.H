/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

Class
    Foam::edge

Description
    An edge is a list of two point labels. The functionality it provides
    supports the discretisation on a 2-D flat mesh.

SourceFiles
    edgeI.H

\*---------------------------------------------------------------------------*/

#ifndef edge_H
#define edge_H

#include "FixedList.H"
#include "pointField.H"
#include "linePointRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class edge;
inline bool operator==(const edge& a, const edge& b);
inline bool operator!=(const edge& a, const edge& b);


/*---------------------------------------------------------------------------*\
                           Class edge Declaration
\*---------------------------------------------------------------------------*/

class edge
:
    public FixedList<label, 2>
{

public:

    // Static data members

        static const char* const typeName;


    // Constructors

        //- Construct null with invalid point labels (-1)
        inline edge();

        //- Construct from components
        inline edge(const label from, const label to);

        //- Construct, optionally sorted with start less-than end
        inline edge(const label from, const label to, const bool doSort);

        //- Construct from FixedList
        inline edge(const FixedList<label, 2>& lst);

        //- Construct, optionally sorted with start less-than end
        inline edge(const FixedList<label, 2>& lst, const bool doSort);

        //- Construct from Istream
        inline edge(Istream& is);


    // Member Functions

        //- Return start vertex label
        inline label start() const;

        //- Return start vertex label
        inline label& start();

        //- Return end vertex label
        inline label end() const;

        //- Return end vertex label
        inline label& end();

        //- Do the edges share a common vertex index?
        inline bool connects(const edge& other) const;

        //- Return vertex common with otherEdge or -1 on failure
        inline label commonVertex(const edge& other) const;

        //- Given one vertex index, return the other one.
        inline label otherVertex(const label index) const;

        //- 'Collapse' edge by marking duplicate point labels.
        //  Duplicates point labels are marked with '-1'.
        //  (the lower vertex is retained).
        //  Return the collapsed size.
        inline label collapse();

        //- Return true if point label is found in edge
        //  No special treatment for '-1'.
        inline bool found(const label index) const;

        //- Return the number of unique, valid (non -1) point labels.
        //  Similar to a HashTable::size().
        inline label count() const;

        //- Insert the index if it did not previously exist on the edge.
        //  Returns true on success. A negative label never inserts.
        //  Similar to a HashTable::insert().
        inline bool insert(const label index);

        //- Remove an existing index from the edge and set its location to '-1'.
        //  Returns true on success. A negative label never removes.
        //  Similar to a HashTable::erase().
        inline bool erase(const label index);


        //- True if the edge is sorted such that start is less-than end
        inline bool sorted() const;

        //- Sort start/end that start is less-than end
        inline void sort();

        //- Flip the edge in-place.
        inline void flip();

        //- Return reverse edge
        inline edge reverseEdge() const;

        //- Return centre (centroid)
        inline point centre(const UList<point>& pts) const;

        //- Return the vector (end - start)
        inline vector vec(const UList<point>& pts) const;

        //- Return the unit vector (end - start)
        inline vector unitVec(const UList<point>& pts) const;

        //- Return scalar magnitude
        inline scalar mag(const UList<point>& pts) const;

        //- Return edge line
        inline linePointRef line(const UList<point>& pts) const;

        //- Compare edges
        //  Returns:
        //  -  0: different
        //  - +1: identical
        //  - -1: same edge, but different orientation
        static inline int compare(const edge& a, const edge& b);


    // Friend Operators

        friend bool operator==(const edge& a, const edge& b);
        friend bool operator!=(const edge& a, const edge& b);
};


//- Hash specialization for hashing edges - a commutative hash value.
//  Hash incrementally.
template<>
inline unsigned Hash<edge>::operator()(const edge& e, unsigned seed) const
{
    unsigned val = seed;

    if (e[0] < e[1])
    {
        val = Hash<label>()(e[0], val);
        val = Hash<label>()(e[1], val);
    }
    else
    {
        val = Hash<label>()(e[1], val);
        val = Hash<label>()(e[0], val);
    }

    return val;
}


//- Hash specialization for hashing edges - a commutative hash value.
//  Hash incrementally.
template<>
inline unsigned Hash<edge>::operator()(const edge& e) const
{
    return Hash<edge>()(e, 0);
}


template<>
inline bool contiguous<edge>()  {return true;}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "edgeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
