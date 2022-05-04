/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::wallFunctionBlenders

Description
    Class to describe the wall-function coefficients being
    used in the wall function boundary conditions.

SourceFiles
    wallFunctionBlenders.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_wallFunctionBlenders_H
#define Foam_wallFunctionBlenders_H

#include "scalarFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;
class Ostream;
template<class EnumType> class Enum;

/*---------------------------------------------------------------------------*\
                    Class wallFunctionBlenders Declaration
\*---------------------------------------------------------------------------*/

class wallFunctionBlenders
{
protected:

    // Protected Enumerations

        //- Options for the blending treatment of viscous and inertial sublayers
        enum blenderType
        {
            STEPWISE,       //!< "Stepwise switch (discontinuous)"
            MAX,            //!< "Maximum value switch (discontinuous)"
            BINOMIAL,       //!< "Binomial blending (smooth)"
            EXPONENTIAL,    //!< "Exponential blending (smooth)"
            TANH            //!< "Tanh blending (smooth)"
        };

        //- Names for blenderType
        static const Enum<blenderType> blenderTypeNames;


    // Protected Data

        //- Blending treatment
        enum blenderType blender_;

        //- Binomial blending exponent being used when
        //- blenderType is blenderType::BINOMIAL
        scalar n_;


public:

    // Constructors

        //- Default construct with default coefficients
        wallFunctionBlenders();

        //- Construct from dictionary and default parameters
        wallFunctionBlenders
        (
            const dictionary& dict,
            const blenderType blender,
            const scalar n
        );


    // Member Functions

    // I-O

        //- Write wall-function blending data as dictionary entries
        void writeEntries(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
