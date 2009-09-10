/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::collisionRecord

Description

    Record of a collision between the particle holding the record and
    the particle with the stored id.

    The access status of the record is to be coded in the
    origProcOfOther data member.  The actual processor is offset by
    +1.  A negative value means that the record has not been accessed,
    positive means that it has.

SourceFiles
    collisionRecordI.H
    collisionRecord.C
    collisionRecordIO.C

\*---------------------------------------------------------------------------*/

#ifndef collisionRecord_H
#define collisionRecord_H

#include "label.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class collisionRecord;
Istream& operator>>(Istream&, collisionRecord&);
Ostream& operator<<(Ostream&, const collisionRecord&);


/*---------------------------------------------------------------------------*\
                         Class collisionRecord Declaration
\*---------------------------------------------------------------------------*/

class collisionRecord
{
    // Private data

        //- Originating processor id of other collision partner
        label origProcOfOther_;

        //- Local particle id on originating processor of other
        //  collision partner
        label origIdOfOther_;

        //- Collision data, stored as if the storing particle was the
        //  first particle (particle A) in the collision.
        vector data_;


public:

    // Constructors

        //- Construct null
        collisionRecord();

        //- Construct from components
        collisionRecord
        (
            label origProcOfOther,
            label origIdOfOther,
            const vector& data = vector::zero
        );

        //- Construct from Istream
        collisionRecord(Istream&);

        //- Construct as copy
        collisionRecord(const collisionRecord&);


    //- Destructor
    ~collisionRecord();


    // Member Functions

        // Access

            //- Return the origProcOfOther data
            inline label origProcOfOther() const;

            //- Return the origIdOfOther data
            inline label origIdOfOther() const;

            //- Return access to the collision data
            inline const vector& collisionData() const;

            //- Return access to the collision data
            inline vector& collisionData();


        // Check

            inline bool match
            (
                label queryOrigProcOfOther,
                label queryOrigIdOfOther
            ) const;

            //- Return the accessed status of the record
            inline bool accessed() const;


        // Edit

            //- Set the accessed property of the record to accessed
            inline void setAccessed();

            //- Set the accessed property of the record to unaccessed
            inline void setUnaccessed();


    // Member Operators

        void operator=(const collisionRecord&);


    // Friend Operators

        friend bool operator==(const collisionRecord&, const collisionRecord&);
        friend bool operator!=(const collisionRecord&, const collisionRecord&);


    // IOstream Operators

        friend Istream& operator>>(Istream&, collisionRecord&);
        friend Ostream& operator<<(Ostream&, const collisionRecord&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "collisionRecordI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
