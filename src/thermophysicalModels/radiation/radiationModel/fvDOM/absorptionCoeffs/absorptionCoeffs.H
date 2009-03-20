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
    Foam::janafThermo

Description
    Absorption coefficients class used in greyMeanAbsorptionEmission and
    wideBandAbsorptionEmission

SourceFiles
    absorptionCoeffs.C

\*---------------------------------------------------------------------------*/

#ifndef absorptionCoeffs_H
#define absorptionCoeffs_H

#include "List.H"
#include "IOstreams.H"
#include "IOdictionary.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                           Class absorptionCoeffs Declaration
\*---------------------------------------------------------------------------*/

class absorptionCoeffs
{

public:

        static const int nCoeffs_ = 6;
        typedef FixedList<scalar, nCoeffs_> coeffArray;

private:

    // Private data

        // Temperature limits of applicability of functions
        scalar Tcommon_;

        scalar Tlow_;

        scalar Thigh_;

        // Polynomio using inverse temperatures
        bool invTemp_;

        coeffArray highACoeffs_;
        coeffArray lowACoeffs_;


    // Private member functions

        //- Check given temperature is within the range of the fitted coeffs
        void checkT(const scalar T) const;


public:

    // Constructors

        //- Construct from Istream
        absorptionCoeffs(Istream&);

        // Null constructor
        absorptionCoeffs()
        {}

    // Destructor

        ~absorptionCoeffs();

    // member functions

        //- Return the coefficients corresponding to the given temperature
        const coeffArray& coeffs(const scalar T) const;

        // Init from a dictionary
        void init(const dictionary&);

        // Init from an Istram
        inline void init(Istream&)
        {
             absorptionCoeffs(Istream);
        }


    // Acces Functions
        inline bool invTemp() const
        {
            return  invTemp_;
        }

        inline scalar Tcommon() const
        {
            return  Tcommon_;
        }

        inline scalar Tlow() const
        {
            return  Tlow_;
        }

        inline scalar Thigh() const
        {
            return  Thigh_;
        }

        inline const coeffArray& highACoeffs() const
        {
            return  highACoeffs_;
        }

        inline const coeffArray& lowACoeffs() const
        {
            return  lowACoeffs_;
        }

};

} // End namespace Foam

} // End namespace radiation

#endif

