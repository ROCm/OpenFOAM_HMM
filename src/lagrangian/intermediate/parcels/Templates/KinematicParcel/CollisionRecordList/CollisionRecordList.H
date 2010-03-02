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
    Foam::CollisionRecordList

Description

SourceFiles
    CollisionRecordListI.H
    CollisionRecordList.C
    CollisionRecordListIO.C

\*---------------------------------------------------------------------------*/

#ifndef CollisionRecordList_H
#define CollisionRecordList_H

#include "DynamicList.H"
#include "CollisionRecord.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class CollisionRecordList Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CollisionRecordList
:
    public DynamicList<CollisionRecord<Type> >
{

public:

    // Constructors

        //- Construct null
        CollisionRecordList();

        //- Construct from Istream
        CollisionRecordList(Istream&);


    //- Destructor
    ~CollisionRecordList();


    // Member Functions

        //- Enquires if the proc and id pair of the other particle are
        //  present in the records.  If so, return access to the
        //  collisionData and mark the CollisionRecord as accessed this
        //  step, if not, create the record and return access to it.
        CollisionRecord<Type>& matchRecord
        (
            label origProcOfOther,
            label origIdOfOther
        );

        //- Update the collision records, deleting any records not
        //  marked as having been accessed, then mark all records as
        //  not accessed ready for the next evaluation
        void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CollisionRecordList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
