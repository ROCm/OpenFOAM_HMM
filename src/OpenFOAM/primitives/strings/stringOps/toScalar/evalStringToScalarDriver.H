/*---------------------------------------------------------------------------*\
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
    Foam::parsing::evalStringToScalar::parseDriver

Description
    Driver for stringOps::toScalar parsing

SourceFiles
    evalStringToScalarDriver.C

\*---------------------------------------------------------------------------*/

#ifndef evalStringToScalarDriver_H
#define evalStringToScalarDriver_H

#include "scalar.H"
#include "className.H"
#include "genericRagelLemonDriver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace parsing
{
namespace evalStringToScalar
{

/*---------------------------------------------------------------------------*\
                         Class parseDriver Declaration
\*---------------------------------------------------------------------------*/

class parseDriver
:
    public genericRagelLemonDriver
{
protected:

    // Protected Data

        //- The result
        scalar value_;


    // Protected Member Functions

        // No copy copy construct
        parseDriver(const parseDriver&) = delete;

        // No copy assignment
        void operator=(const parseDriver&) = delete;


public:

    ClassName("evalStringToScalar::driver");

    // Constructors

        //- Construct null
        parseDriver();


    //- Destructor
    virtual ~parseDriver() = default;


    // Member Functions

        //- Perform parsing on (sub) string
        //  \return evaluated value
        scalar execute
        (
            const std::string& s,
            size_t strPos = 0,
            size_t strLen = std::string::npos
        );

        //- Get value
        scalar value() const
        {
            return value_;
        }

        //- Set value
        void setValue(scalar val)
        {
            value_ = val;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace evalStringToScalar
} // End namespace parsing
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
