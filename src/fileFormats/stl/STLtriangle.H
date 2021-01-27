/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2017 OpenCFD Ltd.
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
    Foam::STLtriangle

Description
    A triangle representation for STL files.

SourceFiles
    STLtriangleI.H

\*---------------------------------------------------------------------------*/

#ifndef STLtriangle_H
#define STLtriangle_H

#include "STLpoint.H"
#include "Istream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class STLtriangle;
inline Ostream& operator<<(Ostream& os, const STLtriangle& tri);

/*---------------------------------------------------------------------------*\
                         Class STLtriangle Declaration
\*---------------------------------------------------------------------------*/

class STLtriangle
{
    // Typedefs

    //- Attribute is 16-bit
    typedef unsigned short STLattrib;

    // Private data
    // NB: The order of the members (1 normal, 3 points, 1 attribute) is
    // essential when reading/writing binary content.

        //- The face normal and the three points defining the triangle.
        //  Some programs may write zero or other junk for the normal.
        STLpoint normal_, a_, b_, c_;

        //- The attribute information could be for 'color' or solid id, etc
        STLattrib attrib_;


public:

    // Constructors

        //- Construct null
        inline STLtriangle();

        //- Construct from components
        inline STLtriangle
        (
            const STLpoint& normal,
            const STLpoint& a,
            const STLpoint& b,
            const STLpoint& c,
            unsigned short attrib
        );

        //- Construct from istream (read binary)
        inline STLtriangle(std::istream& is);


    // Member Functions

        // Access

            inline const STLpoint& normal() const;
            inline const STLpoint& a() const;
            inline const STLpoint& b() const;
            inline const STLpoint& c() const;
            inline unsigned short attrib() const;


        // Read

            //- Read from istream (binary)
            inline void read(std::istream& is);


        // Write

            //- Write to ostream (binary)
            inline void write(std::ostream& os) const;

            //- Write to Ostream (ASCII)
            inline Ostream& print(Ostream& os) const;


            //- Write components to Ostream (ASCII)
            inline static void write
            (
                Ostream& os,
                const vector& norm,
                const point& pt0,
                const point& pt1,
                const point& pt2
            );

            //- Write components to Ostream (ASCII), calculating the normal
            inline static void write
            (
                Ostream& os,
                const point& pt0,
                const point& pt1,
                const point& pt2
            );


    // Ostream Operator

        //- Print triangle contents
        friend Ostream& operator<<(Ostream& os, const STLtriangle& tri);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "STLtriangleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
