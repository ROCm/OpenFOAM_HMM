/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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
    Foam::tabulatedWallFunctions::general

Description
    Computes U+ as a function of Reynolds number by inverting table of
    y+ vs U+

    Example dictionary specification:

        tabulatedWallFunction general;

        // Output table info
        tableName       uPlusWallFunctionData; // Output table name
        log10           yes;                // Re interpreted as log10(Rey)
        dx              0.2;                // Interval log10(Rey)
        x0              -3;                 // Minimum log10(Rey)
        xMax            7;                  // Maximum log10(Rey)

        generalCoeffs
        {
            interpType      linear;         // Interpolation method
            log10YPlus      true;           // y+ values defined as log10(y+)
            log10UPlus      true;           // U+ values defined as log10(y+)
            inputTable
            (
                (yPlusValue0    uPlusValue0)
                ...
                (yPlusValueN    uPlusValueN)
            );

        }


SourceFiles
    general.C

\*---------------------------------------------------------------------------*/

#ifndef general_H
#define general_H

#include "tabulatedWallFunction.H"
#include "NamedEnum.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace tabulatedWallFunctions
{

/*---------------------------------------------------------------------------*\
                        Class general Declaration
\*---------------------------------------------------------------------------*/

class general
:
    public tabulatedWallFunction
{
public:

    // Public data types

        //- Enumeration listing available interpolation types
        enum interpolationType
        {
            itLinear
        };

        static const NamedEnum<interpolationType, 1> interpolationTypeNames_;


protected:

    // Protected data

        //- Type of interpolation to apply when inverting the data set
        interpolationType interpType_;

        //- Input y+ values
        List<scalar> yPlus_;

        //- Input U+ values
        List<scalar> uPlus_;

        //- Are y+ values entered as log10(y+)?
        Switch log10YPlus_;

        //- Are U+ values entered as log10(U+)?
        Switch log10UPlus_;


    // Protected member functions

        //- Invert the table
        virtual void invertTable();

        //- Interpolate
        virtual scalar interpolate
        (
            const scalar xi,
            const scalarList& x,
            const scalarList& fx
        ) const;


public:

    //- Run-time type information
    TypeName("general");


    // Constructors
    general(const dictionary& dict, const polyMesh& mesh);

    //- Destructor
    virtual ~general();


    // Member Functions

        // Access

            //- Return y+ as a function of u+
            virtual scalar yPlus(const scalar uPlus) const;

            //- Return Reynolds number as a function of u+
            virtual scalar Re(const scalar uPlus) const;


        // I-O

            //- Write to Ostream
            virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace tabulatedWallFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
