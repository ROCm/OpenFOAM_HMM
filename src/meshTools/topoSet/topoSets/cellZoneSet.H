/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    Foam::cellZoneSet

Description
    Like cellSet but -reads data from cellZone -updates cellZone when writing.

SourceFiles
    cellZone.C

\*---------------------------------------------------------------------------*/

#ifndef cellZoneSet_H
#define cellZoneSet_H

#include "cellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cellZoneSet Declaration
\*---------------------------------------------------------------------------*/

class cellZoneSet
:
    public cellSet
{
    // Private Data

        const polyMesh& mesh_;

        labelList addressing_;

public:

    //- Runtime type information
    TypeName("cellZoneSet");


    // Constructors

        //- Construct from objectRegistry and name
        cellZoneSet
        (
            const polyMesh& mesh,
            const word& name,
            readOption r=MUST_READ,
            writeOption w=NO_WRITE
        );

        //- Construct from initial size for labelHashSet
        cellZoneSet
        (
            const polyMesh& mesh,
            const word& name,
            const label size,
            writeOption w=NO_WRITE
        );

        //- Copy construct from existing set
        cellZoneSet
        (
            const polyMesh& mesh,
            const word& name,
            const topoSet& set,
            writeOption w=NO_WRITE
        );



    //- Destructor
    virtual ~cellZoneSet() = default;


    // Member functions

        const labelList& addressing() const
        {
            return addressing_;
        }

        labelList& addressing()
        {
            return addressing_;
        }

        //- Sort addressing and make cellSet part consistent with addressing
        void updateSet();

        //- Invert contents.
        //  Insert all members [0,maxLen) which were not in set
        virtual void invert(const label maxLen);

        //- Subset contents. Only elements present in both sets remain.
        virtual void subset(const topoSet& set);

        //- Add elements present in set.
        virtual void addSet(const topoSet& set);

        //- Subtract elements present in set.
        virtual void subtractSet(const topoSet& set);

        //- Sync cellSet across coupled patches; update cellZone from cellSet
        virtual void sync(const polyMesh& mesh);

        //- Write maxLen items with label and coordinates.
        virtual void writeDebug
        (
            Ostream& os,
            const primitiveMesh&,
            const label maxLen
        ) const;

        //- Write cellZone using stream options
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool valid
        ) const;

        //- Update any stored data for new labels
        virtual void updateMesh(const mapPolyMesh& morphMap);

        //- Return max index+1.
        virtual label maxSize(const polyMesh& mesh) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
