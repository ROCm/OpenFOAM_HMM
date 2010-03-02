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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline bool Foam::CollisionRecord<Type>::match
(
    label queryOrigProcOfOther,
    label queryOrigIdOfOther
) const
{
    return
    (
        queryOrigProcOfOther == origProcOfOther()
     && queryOrigIdOfOther == origIdOfOther()
    );
}


template<class Type>
inline Foam::label Foam::CollisionRecord<Type>::origProcOfOther() const
{
    return mag(origProcOfOther_) - 1;
}


template<class Type>
inline Foam::label Foam::CollisionRecord<Type>::origIdOfOther() const
{
    return origIdOfOther_;
}


template<class Type>
inline const Type&
Foam::CollisionRecord<Type>::collisionData() const
{
    return data_;
}


template<class Type>
inline Type& Foam::CollisionRecord<Type>::collisionData()
{
    return data_;
}


template<class Type>
inline bool Foam::CollisionRecord<Type>::accessed() const
{
    return pos(origProcOfOther_);
}


template<class Type>
inline void Foam::CollisionRecord<Type>::setAccessed()
{
    origProcOfOther_ = origProcOfOther() + 1;
}


template<class Type>
inline void Foam::CollisionRecord<Type>::setUnaccessed()
{
    origProcOfOther_ = -(origProcOfOther() + 1);
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

template<class Type>
inline bool Foam::operator==
(
    const CollisionRecord<Type>& a,
    const CollisionRecord<Type>& b
)
{
    return
    (
        a.origProcOfOther_ == b.origProcOfOther_
     && a.origIdOfOther_ == b.origIdOfOther_
     && a.data_ == b.data_
    );
}

template<class Type>
inline bool Foam::operator!=
(
    const CollisionRecord<Type>& a,
    const CollisionRecord<Type>& b
)
{
    return !(a == b);
}


// ************************************************************************* //
