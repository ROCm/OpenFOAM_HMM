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
    Foam::CollisionRecord

Description

    Record of a collision between the particle holding the record and
    the particle with the stored id.

    The access status of the record is to be coded in the
    origProcOfOther data member.  The actual processor is offset by
    +1.  A negative value means that the record has not been accessed,
    positive means that it has.

SourceFiles
    CollisionRecordI.H
    CollisionRecord.C
    CollisionRecordIO.C

\*---------------------------------------------------------------------------*/

#ifndef CollisionRecord_H
#define CollisionRecord_H

#include "label.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
template<class Type>
class CollisionRecord;

template<class Type>
inline bool operator==
(
    const CollisionRecord<Type>&,
    const CollisionRecord<Type>&
);

template<class Type>
inline bool operator!=
(
    const CollisionRecord<Type>&,
    const CollisionRecord<Type>&
);

template<class Type>
Istream& operator>>(Istream&, CollisionRecord<Type>&);

template<class Type>
Ostream& operator<<(Ostream&, const CollisionRecord<Type>&);


/*---------------------------------------------------------------------------*\
                         Class CollisionRecord Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CollisionRecord
{
    // Private data

        //- Originating processor id of other collision partner
        label origProcOfOther_;

        //- Local particle id on originating processor of other
        //  collision partner
        label origIdOfOther_;

        //- Collision data, stored as if the storing particle was the
        //  first particle (particle A) in the collision.
        Type data_;


public:

    // Constructors

        //- Construct null
        CollisionRecord();

        //- Construct from components
        CollisionRecord
        (
            label origProcOfOther,
            label origIdOfOther,
            const Type& data = pTraits<Type>::zero
        );

        //- Construct from Istream
        CollisionRecord(Istream&);

        //- Construct as copy
        CollisionRecord(const CollisionRecord&);


    //- Destructor
    ~CollisionRecord();


    // Member Functions


        // Access

            //- Return the origProcOfOther data
            inline label origProcOfOther() const;

            //- Return the origIdOfOther data
            inline label origIdOfOther() const;

            //- Return access to the collision data
            inline const Type& collisionData() const;

            //- Return access to the collision data
            inline Type& collisionData();


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

        void operator=(const CollisionRecord&);


    // Friend Operators

        friend bool operator== <Type>
        (
            const CollisionRecord<Type>&,
            const CollisionRecord<Type>&
        );

        friend bool operator!= <Type>
        (
            const CollisionRecord<Type>&,
            const CollisionRecord<Type>&
        );


    // IOstream Operators

        friend Istream& operator>> <Type>
        (
            Istream&,
            CollisionRecord<Type>&
        );

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const CollisionRecord<Type>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CollisionRecordI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CollisionRecord.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
