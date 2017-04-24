/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::orientedType

Description
    Class to determine the 'oriented' status of surface fields

SourceFiles
    orientedType.C

\*---------------------------------------------------------------------------*/

#ifndef orientedType_H
#define orientedType_H

#include "Istream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class orientedType;

Istream& operator>>(Istream&, orientedType&);

Ostream& operator<<(Ostream&, const orientedType&);

/*---------------------------------------------------------------------------*\
                           Class dimensioned Declaration
\*---------------------------------------------------------------------------*/

class orientedType
{
    // Private data

        //- Oriented flag
        bool oriented_;


public:


    // Constructors

        //- Null constructor - flag initialised to false
        orientedType();

        //- Copy constructor
        orientedType(const orientedType& of);

        //- Construct from bool
        orientedType(const bool oriented);

        //- Construct from Istream
        orientedType(Istream& is);


    // Member functions

        //- Return non-const reference to the oriented flag
        bool& oriented();

        //- Return const reference to the oriented flag
        bool oriented() const;


    // Member operators

        void operator=(const orientedType& of);

        void operator+=(const orientedType& of);
        void operator-=(const orientedType& of);
        void operator*=(const orientedType& of);
        void operator/=(const orientedType& of);
        void operator*=(const scalar s);
        void operator/=(const scalar s);
        bool operator()() const;


    // IOstream operators

        friend Istream& operator>>(Istream& is, orientedType& of);

        friend Ostream& operator<<(Ostream& os, const orientedType& of);
};


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

orientedType max(const orientedType& of1, const orientedType& of2);
orientedType min(const orientedType& of1, const orientedType& of2);
orientedType cmptMultiply(const orientedType& of1, const orientedType& of2);
orientedType cmptDivide(const orientedType& of1, const orientedType& of);
orientedType cmptAv(const orientedType& of);


orientedType pow(const orientedType& of, const scalar r);
orientedType sqr(const orientedType& of);
orientedType pow3(const orientedType& of);
orientedType pow4(const orientedType& of);
orientedType pow5(const orientedType& of);
orientedType pow6(const orientedType& of);
orientedType pow025(const orientedType& of);


orientedType sqrt(const orientedType& of);
orientedType cbrt(const orientedType& of);
orientedType magSqr(const orientedType& of);
orientedType mag(const orientedType& of);
orientedType sign(const orientedType& of);
orientedType pos(const orientedType& of);
orientedType neg(const orientedType& of);
orientedType posPart(const orientedType& of);
orientedType negPart(const orientedType& of);
orientedType inv(const orientedType& of);


orientedType trans(const orientedType& of);
orientedType atan2(const orientedType& of1, const orientedType& of2);
orientedType transform(const orientedType& of);

orientedType operator-(const orientedType& of);
orientedType operator*(const scalar s, const orientedType& of);
orientedType operator/(const orientedType& of, const scalar s);

orientedType operator+(const orientedType& of1, const orientedType& of2);
orientedType operator-(const orientedType& of1, const orientedType& of2);
orientedType operator/(const orientedType& of1, const orientedType& of2);
orientedType operator*(const orientedType& of1, const orientedType& of2);
orientedType operator^(const orientedType& of1, const orientedType& of2);
orientedType operator&(const orientedType& of1, const orientedType& of2);
orientedType operator&&(const orientedType& of1, const orientedType& of2);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
