/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class orientedType;

Istream& operator>>(Istream& is, orientedType& ot);
Ostream& operator<<(Ostream& os, const orientedType& ot);

/*---------------------------------------------------------------------------*\
                        Class orientedType Declaration
\*---------------------------------------------------------------------------*/

class orientedType
{
public:

    // Public Data Types

        //- Enumeration defining oriented flags
        enum orientedOption
        {
            UNKNOWN = 0,
            ORIENTED = 1,
            UNORIENTED = 2
        };

        //- Named enumerations for oriented flags
        static const Enum<orientedOption> orientedOptionNames;


private:

    // Private Data

        //- Oriented type
        orientedOption oriented_;


public:

    // Constructors

        //- Default construct as "UNKNOWN"
        orientedType();

        //- Copy construct
        orientedType(const orientedType& ot);

        //- Construct from bool
        explicit orientedType(const bool oriented);

        //- Construct from Istream
        explicit orientedType(Istream& is);


    // Member Functions

        //- Return true if can operate on this pair of oriented types
        static bool checkType
        (
            const orientedType& ot1,
            const orientedType& ot2
        );

        //- Return non-const reference to the oriented flag
        orientedOption& oriented() noexcept;

        //- Return the oriented flag
        orientedOption oriented() const noexcept;

        //- Set the oriented flag
        void setOriented(const bool oriented = true) noexcept;

        //- Read the "oriented" state from dictionary
        void read(const dictionary& dict);

        //- Write the "oriented" flag entry (if ORIENTED)
        //  \return True if entry was written
        bool writeEntry(Ostream& os) const;


    // Member Operators

        void operator=(const orientedType& ot);

        void operator+=(const orientedType& ot);
        void operator-=(const orientedType& ot);
        void operator*=(const orientedType& ot);
        void operator/=(const orientedType& ot);
        void operator*=(const scalar s);
        void operator/=(const scalar s);
        bool operator()() const;


    // IOstream Operators

        friend Istream& operator>>(Istream& is, orientedType& ot);
        friend Ostream& operator<<(Ostream& os, const orientedType& ot);
};


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

orientedType max(const orientedType& ot1, const orientedType& ot2);
orientedType min(const orientedType& ot1, const orientedType& ot2);
orientedType cmptMultiply(const orientedType& ot1, const orientedType& ot2);
orientedType cmptDivide(const orientedType& ot1, const orientedType& ot);
orientedType cmptAv(const orientedType& ot);


orientedType pow(const orientedType& ot, const scalar r);
orientedType sqr(const orientedType& ot);
orientedType pow3(const orientedType& ot);
orientedType pow4(const orientedType& ot);
orientedType pow5(const orientedType& ot);
orientedType pow6(const orientedType& ot);
orientedType pow025(const orientedType& ot);


orientedType sqrt(const orientedType& ot);
orientedType cbrt(const orientedType& ot);
orientedType magSqr(const orientedType& ot);
orientedType mag(const orientedType& ot);
orientedType sign(const orientedType& ot);
orientedType pos(const orientedType& ot);
orientedType pos0(const orientedType& ot);
orientedType neg(const orientedType& ot);
orientedType neg0(const orientedType& ot);
orientedType posPart(const orientedType& ot);
orientedType negPart(const orientedType& ot);
orientedType inv(const orientedType& ot);


orientedType trans(const orientedType& ot);
orientedType atan2(const orientedType& ot1, const orientedType& ot2);
orientedType hypot(const orientedType& ot1, const orientedType& ot2);
orientedType transform(const orientedType& ot);

orientedType operator-(const orientedType& ot);
orientedType operator*(const scalar s, const orientedType& ot);
orientedType operator/(const orientedType& ot, const scalar s);

orientedType operator+(const orientedType& ot1, const orientedType& ot2);
orientedType operator-(const orientedType& ot1, const orientedType& ot2);
orientedType operator/(const orientedType& ot1, const orientedType& ot2);
orientedType operator*(const orientedType& ot1, const orientedType& ot2);
orientedType operator^(const orientedType& ot1, const orientedType& ot2);
orientedType operator&(const orientedType& ot1, const orientedType& ot2);
orientedType operator&&(const orientedType& ot1, const orientedType& ot2);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
