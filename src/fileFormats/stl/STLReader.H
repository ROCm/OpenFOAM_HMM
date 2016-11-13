/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::fileFormats::STLReader

Description
    Internal class used by the STLsurfaceFormat

SourceFiles
    STLReader.C
    STLsurfaceFormatASCII.L

\*---------------------------------------------------------------------------*/

#ifndef STLReader_H
#define STLReader_H

#include "STLCore.H"
#include "labelledTri.H"
#include "IFstream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                   Class fileFormats::STLReader Declaration
\*---------------------------------------------------------------------------*/

class STLReader
:
    public STLCore
{
    // Private Data

        bool sorted_;

        //- The points supporting the facets
        pointField points_;

        //- The zones associated with the faces
        List<label> zoneIds_;

        //- The solid names, in the order of their first appearance
        List<word> names_;

        //- The solid count, in the order of their first appearance
        List<label> sizes_;


    // Private Member Functions

        //- Read ASCII
        bool readASCII(const fileName&);

        //- Read BINARY
        bool readBINARY(const fileName&);

        //- Read ASCII or BINARY
        bool readFile(const fileName&, const STLFormat&);


        //- Disallow default bitwise copy construct
        STLReader(const STLReader&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const STLReader&) = delete;


public:

    // Constructors

        //- Read from file, filling in the information
        STLReader(const fileName&);

        //- Read from file, filling in the information.
        //  Manually selected choice of ascii/binary/detect.
        STLReader(const fileName&, const STLFormat&);


    //- Destructor
    ~STLReader();


    // Member Functions

        //- File read was already sorted
        bool sorted() const
        {
            return sorted_;
        }

        //- Flush all values
        void clear()
        {
            sorted_ = true;
            points_.clear();
            zoneIds_.clear();
            names_.clear();
            sizes_.clear();
        }

        //- Return full access to the points
        pointField& points()
        {
            return points_;
        }

        //- Return full access to the zoneIds
        List<label>& zoneIds()
        {
            return zoneIds_;
        }

        //- The list of solid names in the order of their first appearance
        List<word>& names()
        {
            return names_;
        }

        //- The list of solid sizes in the order of their first appearance
        List<label>& sizes()
        {
            return sizes_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
