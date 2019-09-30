/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::parsing::evalStringToScalar::parser

Description
    Interface to lemon parser to simple string to scalar evaluation, which
    is used by stringOps::toScalar

\*---------------------------------------------------------------------------*/

#ifndef evalStringToScalarParser_H
#define evalStringToScalarParser_H

namespace Foam
{
namespace parsing
{
namespace evalStringToScalar
{

// Forward Declarations
class parseDriver;


/*---------------------------------------------------------------------------*\
                           Class parser Declaration
\*---------------------------------------------------------------------------*/

class parser
{
    // Private Data

        //- The lemon parser (demand-driven)
        void* lemon_;

public:

    // Constructors

        //- Construct null
        parser() : lemon_(nullptr) {}


    //- Destructor, delete lemon parser
    ~parser()
    {
        stop();
    }


    // Member Functions

        //- Start parsing, with given driver context
        void start(parseDriver& driver);

        //- Stop parsing, freeing the allocated parser
        void stop();

        //- Push token/value to parser
        void parse(int tokenId, Foam::scalar val);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace evalStringToScalar
} // End namespace parsing
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
