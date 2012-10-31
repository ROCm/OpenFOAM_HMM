/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Class
    Foam::solidSpecie

Description
    Base class of the solid thermophysical property types.

SourceFiles
    solidSpecieI.H
    solidSpecie.C

\*---------------------------------------------------------------------------*/

#ifndef solidSpecie_H
#define solidSpecie_H

#include "word.H"
#include "scalar.H"
#include "dictionary.H"
#include "specie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

/*---------------------------------------------------------------------------*\
                           Class solidSpecie Declaration
\*---------------------------------------------------------------------------*/

class solidSpecie
:
    public specie
{

    // Private Member Functions

        //- Construct from components without name
        inline solidSpecie(const scalar nMoles, const scalar molWeight);


public:

    //- Runtime type information
    ClassName("solidSpecie");


    // Constructors

        //- Construct from components with name
        inline solidSpecie
        (
            const word& name,
            const scalar nMoles,
            const scalar molWeight
        );

        //- Construct as copy
        inline solidSpecie(const solidSpecie&);

        //- Construct as named copy
        inline solidSpecie(const word& name, const solidSpecie&);

        //- Construct from dictionary
        solidSpecie(const dictionary& dict);


    // Member Functions


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator=(const solidSpecie&);

        inline void operator+=(const solidSpecie&);
        inline void operator-=(const solidSpecie&);

        inline void operator*=(const scalar);


    // Friend operators

        inline friend solidSpecie operator+
        (
            const solidSpecie&, const solidSpecie&
        );
        inline friend solidSpecie operator-
        (
            const solidSpecie&, const solidSpecie&
        );

        inline friend solidSpecie operator*(const scalar, const solidSpecie&);

        inline friend solidSpecie operator==
        (
            const solidSpecie&, const solidSpecie&
        );

    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const solidSpecie&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidSpecieI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
