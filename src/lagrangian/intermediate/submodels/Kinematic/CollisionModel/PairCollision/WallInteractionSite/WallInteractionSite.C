/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "WallInteractionSite.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::WallInteractionSite<Type>::WallInteractionSite()
:
    patchI_(),
    wallData_()
{}


template<class Type>
Foam::WallInteractionSite<Type>::WallInteractionSite
(
    label patchI,
    const Type& wallData
)
:
    patchI_(patchI),
    wallData_(wallData)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::WallInteractionSite<Type>::~WallInteractionSite()
{}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type>
bool Foam::WallInteractionSite<Type>::operator==
(
    const WallInteractionSite<Type>& rhs
) const
{
    return patchI_ == rhs.patch_ && wallData_ == rhs.wallData_;
}


template<class Type>
bool Foam::WallInteractionSite<Type>::operator!=
(
    const WallInteractionSite<Type>& rhs
) const
{
    return !(*this == rhs);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Istream& Foam::operator>>
(
    Istream& is,
    WallInteractionSite<Type>& wIS
)
{
    is  >> wIS.patchI_ >> wIS.wallData_;

    // Check state of Istream
    is.check
    (
        "Foam::Istream& Foam::operator>>"
        "(Foam::Istream&, Foam::WallInteractionSite<Type>&)"
    );

    return is;
}


template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const WallInteractionSite<Type>& wIS
)
{
    os  << wIS.patchI_ << token::SPACE << wIS.wallData_;

    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::operator<<"
        "(Ostream&, const WallInteractionSite<Type>&)"
    );

    return os;
}


// ************************************************************************* //
