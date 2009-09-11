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

\*---------------------------------------------------------------------------*/

#include "CollisionRecordList.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::CollisionRecordList<Type>::CollisionRecordList()
:
    DynamicList<CollisionRecord<Type> >()
{}


template<class Type>
Foam::CollisionRecordList<Type>::CollisionRecordList(Istream& is)
:
    DynamicList<CollisionRecord<Type> >(is)
{
    // Check state of Istream
    is.check
    (
        "Foam::CollisionRecordList<Type>::CollisionRecordList(Foam::Istream&)"
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * /

template<class Type>
Foam::CollisionRecordList<Type>::~CollisionRecordList()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Foam::CollisionRecord<Type>& Foam::CollisionRecordList<Type>::matchRecord
(
    label origProcOfOther,
    label origIdOfOther
)
{
    // Returning the first record that matches the particle
    // identifiers.  Two records with the same identification is not
    // supported.

    forAll(*this, i)
    {
        CollisionRecord<Type>& cR = (*this)[i];

        if (cR.match(origProcOfOther, origIdOfOther))
        {
            cR.setAccessed();

            return cR;
        }
    }

    // Record not found, create a new one and return it as the last
    // member of the list.  The status of the record will be accessed
    // by construction.

    append(CollisionRecord<Type>(origProcOfOther, origIdOfOther));

    return (*this)[this->size() - 1];
}


template<class Type>
void Foam::CollisionRecordList<Type>::update()
{
    DynamicList<CollisionRecord<Type> > updatedRecords;

    forAll(*this, i)
    {
        if ((*this)[i].accessed())
        {
            (*this)[i].setUnaccessed();

            updatedRecords.append((*this)[i]);
        }
    }

    DynamicList<CollisionRecord<Type> >::operator=(updatedRecords);
}

// ************************************************************************* //
