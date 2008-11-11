/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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
    Foam::fileFormats::FTRfileFormat

Description
    Reading of Foam Trisurface Format

SourceFiles
    FTRfileFormat.C

\*---------------------------------------------------------------------------*/

#ifndef FTRfileFormat_H
#define FTRfileFormat_H

#include "Ostream.H"
#include "OFstream.H"
#include "keyedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                       Class FTRfileFormat Declaration
\*---------------------------------------------------------------------------*/

class FTRfileFormat
:
    public keyedSurface
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        FTRfileFormat(const FTRfileFormat&);

        //- Disallow default bitwise assignment
        void operator=(const FTRfileFormat&);


    //- read compatibility for ftr format
    class ftrPatch
    {
        //- Name of patch
        word name_;

        //- Type of patch
        word type_;

    public:

        const word& name() const
        {
            return name_;
        }

        friend Istream& operator>>(Istream& is, ftrPatch& p)
        {
            is >> p.name_ >> p.type_;
            return is;
        }
    };

public:

    // Constructors

        //- Construct null
        FTRfileFormat();

        //- Construct from file name
        FTRfileFormat(const fileName&, const bool triangulate=false);

    // Selectors

        //- Read file and return keyedSurface
        static autoPtr<keyedSurface> New
        (
            const fileName& fName,
            const bool triangulate=false
        )
        {
            return autoPtr<keyedSurface>
            (
                new FTRfileFormat(fName,true)
            );
        }

    // Member Functions

        //- Read from file
        virtual bool read(const fileName&, const bool triangulate);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
