/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::instant

Description
    An instant of time. Contains the time value and name.

SourceFiles
    instant.C

\*---------------------------------------------------------------------------*/

#ifndef instant_H
#define instant_H

#include "word.H"
#include "Instant.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class instant Declaration
\*---------------------------------------------------------------------------*/

class instant
:
    public Instant<word>
{
public:

    // Static Data Members

        static const char* const typeName;


    // Constructors

        //- Copy and move construct from components
        using Instant<word>::Instant;

        //- Construct null, with time-value = 0.
        instant() = default;

        //- Copy construct
        instant(const instant&) = default;

        //- Move construct
        instant(instant&&) = default;

        //- Construct from timeValue, auto generating the name
        explicit instant(scalar timeValue);

        //- Construct from timeName, parsing timeName for a value
        explicit instant(const word& timeName);

        //- Construct from timeName, parsing timeName for a value
        explicit instant(word&& timeName);


    // Operators

        //- Copy assignment
        instant& operator=(const instant&) = default;

        //- Move assignment
        instant& operator=(instant&&) = default;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
