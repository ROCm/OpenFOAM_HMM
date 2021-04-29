/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::Function1Types::LimitRange

Description
    Function1 wrapper that limits the input range of another Function1

    Example usage for limiting a polynomial:
    \verbatim
        <entryName>
        {
            type limitRange;
            min  0.4;
            max  1.4;

            value polynomial
            (
                (5 1)
                (-2 2)
                (-2 3)
                (1 4)
            );
        }
    \endverbatim

    Here the return value will be:
    - poly(0.4) for x <= 0.4;
    - poly(1.4) for x >= 1.4; and
    - poly(x) for 0.4 < x < 1.4.


    Example usage for limiting a file-based table:
    \verbatim
        <entryName>
        {
            type limitRange;
            min  0.4;
            max  1.4;

            value
            {
                type table;
                file "<system>/fanCurve.txt";
            }
        }
    \endverbatim

    Where:
    \table
        Property | Description                                  | Required
        min      | Minimum input value                          | yes
        max      | Maximum input value                          | yes
        value    | Function of type Function1<Type>             | yes
    \endtable

SourceFiles
    LimitRange.C

\*---------------------------------------------------------------------------*/

#ifndef Function1Types_LimitRange_H
#define Function1Types_LimitRange_H

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{

/*---------------------------------------------------------------------------*\
                          Class LimitRange Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class LimitRange
:
    public Function1<Type>
{
    // Private Data

        //- Minimum input value
        scalar min_;

        //- Maximum input value
        scalar max_;

        //- Value function
        autoPtr<Function1<Type>> value_;


    // Private Member Functions

        //- Read the coefficients from the given dictionary
        void read(const dictionary& coeffs);


public:

    //- Runtime type information
    TypeName("limitRange");


    // Generated Methods

        //- No copy assignment
        void operator=(const LimitRange<Type>&) = delete;


    // Constructors

        //- Construct from entry name and dictionary
        LimitRange(const word& entryName, const dictionary& dict);

        //- Copy construct
        explicit LimitRange(const LimitRange<Type>& rhs);


    //- Destructor
    virtual ~LimitRange() = default;


    // Member Functions

        //- Return value for time t
        virtual inline Type value(const scalar t) const;

        //- Integrate between two (scalar) values
        virtual inline Type integrate(const scalar x1, const scalar x2) const;

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;

        //- Write coefficient entries in dictionary format
        void writeEntries(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "LimitRangeI.H"

#ifdef NoRepository
    #include "LimitRange.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
