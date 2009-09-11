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

#include "CollisionRecord.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::CollisionRecord<Type>::CollisionRecord()
:
    origProcOfOther_(-VGREAT),
    origIdOfOther_(-VGREAT),
    data_(pTraits<Type>::zero)
{}


template<class Type>
Foam::CollisionRecord<Type>::CollisionRecord
(
    label origProcOfOther,
    label origIdOfOther,
    const Type& data
)
:
    origProcOfOther_(origProcOfOther + 1),
    origIdOfOther_(origIdOfOther),
    data_(data)
{}


template<class Type>
Foam::CollisionRecord<Type>::CollisionRecord(const CollisionRecord<Type>& cR)
:
    origProcOfOther_(cR.origProcOfOther() + 1),
    origIdOfOther_(cR.origIdOfOther_),
    data_(cR.data_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::CollisionRecord<Type>::~CollisionRecord()
{}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type>
void Foam::CollisionRecord<Type>::operator=(const CollisionRecord<Type>& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn
        (
            "Foam::CollisionRecord<Type>::operator="
            "(const Foam::CollisionRecord<Type>&)"
        )
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    origProcOfOther_ = rhs.origProcOfOther_;
    origIdOfOther_ = rhs.origIdOfOther_;
    data_ = rhs.data_;
}


// * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * * * //

#include "CollisionRecordIO.C"


// ************************************************************************* //
