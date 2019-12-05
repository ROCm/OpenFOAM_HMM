/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2018 OpenCFD Ltd.
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
    Foam::ZoneMesh

Description
    A list of mesh zones.

SourceFiles
    ZoneMesh.C

\*---------------------------------------------------------------------------*/

#ifndef ZoneMesh_H
#define ZoneMesh_H

#include "List.H"
#include "regIOobject.H"
#include "pointField.H"
#include "Map.H"
#include "bitSet.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations

template<class ZoneType, class MeshType> class ZoneMesh;

template<class ZoneType, class MeshType>
Ostream& operator<<(Ostream& os, const ZoneMesh<ZoneType, MeshType>& zones);


/*---------------------------------------------------------------------------*\
                           Class ZoneMesh Declaration
\*---------------------------------------------------------------------------*/

template<class ZoneType, class MeshType>
class ZoneMesh
:
    public PtrList<ZoneType>,
    public regIOobject
{
    // Private data

        //- Reference to mesh
        const MeshType& mesh_;

        //- Map of zone labels for given element
        mutable Map<label>* zoneMapPtr_;


    // Private Member Functions

        //- Read if IOobject flags set. Return true if read.
        bool read();

        //- Create zone map
        void calcZoneMap() const;

        //- Templated implementation for names()
        template<class UnaryMatchPredicate>
        static wordList namesImpl
        (
            const PtrList<ZoneType>& list,
            const UnaryMatchPredicate& matcher,
            const bool doSort
        );

        //- Templated implementation for indices()
        template<class UnaryMatchPredicate>
        static labelList indicesImpl
        (
            const PtrList<ZoneType>& list,
            const UnaryMatchPredicate& matcher
        );

        //- Templated implementation for findIndex()
        template<class UnaryMatchPredicate>
        static label findIndexImpl
        (
            const PtrList<ZoneType>& list,
            const UnaryMatchPredicate& matcher
        );


        //- No copy construct
        ZoneMesh(const ZoneMesh&) = delete;

        //- No copy assignment
        void operator=(const ZoneMesh<ZoneType, MeshType>&) = delete;


public:

    //- Debug switch to disallow the use of generic zones
    static int disallowGenericZones;


    // Constructors

        //- Read constructor given IOobject and a MeshType reference
        ZoneMesh
        (
            const IOobject& io,
            const MeshType& mesh
        );

        //- Construct given size
        ZoneMesh
        (
            const IOobject& io,
            const MeshType& mesh,
            const label size
        );

        //- Construct given a PtrList
        ZoneMesh
        (
            const IOobject& io,
            const MeshType& mesh,
            const PtrList<ZoneType>& pzm
        );


    //- Destructor
    ~ZoneMesh();


    // Member Functions

        //- Return the mesh reference
        const MeshType& mesh() const
        {
            return mesh_;
        }

        //- Map of zones containing zone index for all zoned elements
        //  Return -1 if the object is not in the zone
        const Map<label>& zoneMap() const;

        //- Given a global object index, return the zone it is in.
        //  If object does not belong to any zones, return -1
        label whichZone(const label objectIndex) const;

        //- Return a list of zone types
        wordList types() const;

        //- A list of the zone names
        wordList names() const;

        //- A list of zone names satisfying the input matcher
        wordList names(const wordRe& matcher) const;

        //- A list of zone names satisfying the input matchers
        wordList names(const wordRes& matcher) const;

        //- Sorted list of the zone names
        wordList sortedNames() const;

        //- Sorted list of zone names satisfying the input matcher
        wordList sortedNames(const wordRe& matcher) const;

        //- Sorted list of zone names satisfying the input matchers
        wordList sortedNames(const wordRes& matcher) const;


        //- Return zone indices for all matches
        labelList indices(const keyType& key) const;

        //- Return zone indices for all matches
        labelList indices(const wordRes& matcher) const;


        //- Return zone index for the first match, return -1 if not found
        label findIndex(const keyType& key) const;

        //- Return zone index for the first match, return -1 if not found
        label findIndex(const wordRes& matcher) const;


        //- Find zone index given a name, return -1 if not found
        label findZoneID(const word& zoneName) const;


        //- Return all elements (cells, faces, points) contained in the
        //- listed zones.
        //  The bitSet is empty (zero-size) if there are no elements matched
        //  anywhere.
        bitSet selection(const labelUList& zoneIds) const;

        //- Return all elements (cells, faces, points) that match the zone
        //- specification as a bitSet.
        //  The bitSet is empty (zero-size) if there are no elements matched
        //  anywhere.
        bitSet selection(const keyType& key) const;

        //- Return all elements (cells, faces, points) that match the zone
        //- specification as a bitSet.
        //  The bitSet is empty (zero-size) if there are no elements matched
        //  anywhere.
        bitSet selection(const wordRes& matcher) const;


        //- Lookup zone by name and return const pointer, nullptr on error.
        const ZoneType* zonePtr(const word& zoneName) const;

        //- Lookup zone by name and return pointer, nullptr on error.
        ZoneType* zonePtr(const word& zoneName);


        //- Clear addressing
        void clearAddressing();

        //- Clear the zones
        void clear();

        //- Check zone definition. Return true if in error.
        bool checkDefinition(const bool report = false) const;

        //- Check whether all procs have all zones and in same order.
        //  \return True if any errors.
        bool checkParallelSync(const bool report = false) const;

        //- Correct zone mesh after moving points
        void movePoints(const pointField& pts);

        //- writeData member function required by regIOobject
        bool writeData(Ostream& os) const;


    // Member Operators

        //- Return const and non-const reference to zone by index.
        using PtrList<ZoneType>::operator[];

        //- Return const reference to zone by name.
        //  Fatal if the zone does not exist.
        const ZoneType& operator[](const word& zoneName) const;

        //- Return reference to an existing zone by name
        //  Fatal if the zone does not exist.
        ZoneType& operator[](const word& zoneName);

        //- Find an existing zone by name or create a new empty one
        //- if required.
        //
        //  To determine if the zone already existed or was newly created,
        //  it will be necessary to add additional logic in the caller.
        //  For example,
        //  \code
        //      const label nOrig = zones.size();
        //
        //      ZoneType& zn = zones("zoneName");
        //
        //      if (nOrig == zones.size()) { existing... } else { new... }
        //  \endcode
        //  \param zoneName the name of the zone
        //  \param verbose report if an existing zone was selected or
        //      a new zone was created.
        //  \return non-const reference to the existing or new zone
        ZoneType& operator()(const word& zoneName, const bool verbose=false);


    // Ostream operator

        friend Ostream& operator<< <ZoneType, MeshType>
        (
            Ostream& os,
            const ZoneMesh<ZoneType, MeshType>& zones
        );


    // Housekeeping

        //- Identical to the indices() method (AUG-2018)
        labelList findIndices(const keyType& key) const
        {
            return indices(key);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ZoneMesh.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
