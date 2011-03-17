/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::infiniteReactionRate::infiniteReactionRate()
{}


inline Foam::infiniteReactionRate::infiniteReactionRate
(
    const speciesTable&,
    const dictionary&
)
{}


inline void Foam::infiniteReactionRate::write(Ostream& os) const
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::infiniteReactionRate::operator()
(
    const scalar,
    const scalar,
    const scalarField&
) const
{
    return (1);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const infiniteReactionRate& rr
)
{
    os  << token::BEGIN_LIST
        << token::END_LIST;
    return os;
}


// ************************************************************************* //
