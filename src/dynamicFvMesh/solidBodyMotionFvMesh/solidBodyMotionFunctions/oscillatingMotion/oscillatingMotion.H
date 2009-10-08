/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

Class
    Foam::solidBodyMotionFunctions::oscillatingMotion

Description
    SolidBodyMotionFvMesh 6DoF motion function. Oscillating displacement.

SourceFiles
    oscillatingMotion.C

\*---------------------------------------------------------------------------*/

#ifndef oscillatingMotion_H
#define oscillatingMotion_H

#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class oscillatingMotion Declaration
\*---------------------------------------------------------------------------*/

class oscillatingMotion
:
    public solidBodyMotionFunction
{
    // Private data

        //- Amplitude
        vector amplitude_;

        //- Radial velocity
        scalar omega_;


    // Private Member Functions

        //- Disallow copy construct
        oscillatingMotion(const oscillatingMotion&);

        //- Disallow default bitwise assignment
        void operator=(const oscillatingMotion&);


public:

    //- Runtime type information
    TypeName("oscillatingMotion");


    // Constructors

        //- Construct from components
        oscillatingMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );


    // Destructor

        virtual ~oscillatingMotion();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
