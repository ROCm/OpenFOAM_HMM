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
    Foam::fileFormats::STLsurfaceFormatCore

Description
    Internal class used by the STLsurfaceFormat

SourceFiles
    STLsurfaceFormatCore.C
    STLsurfaceFormatASCII.L

\*---------------------------------------------------------------------------*/

#ifndef STLsurfaceFormatCore_H
#define STLsurfaceFormatCore_H

#include "STLtriangle.H"
#include "triFace.H"
#include "gzstream.h"
#include "IFstream.H"
#include "Ostream.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                   Class STLsurfaceFormatCore Declaration
\*---------------------------------------------------------------------------*/

class STLsurfaceFormatCore
{
    // Private Data

        bool binary_;

        //- The points supporting the facets
        pointField points_;

        //- The regions associated with the faces
        List<label> regions_;

        label maxRegionId_;

        HashTable<label> groupToPatch_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        STLsurfaceFormatCore(const STLsurfaceFormatCore&);

        //- Disallow default bitwise assignment
        void operator=(const STLsurfaceFormatCore&);

        //- Determine the file type
        static int detectBINARY(const fileName&);

        //- Read ASCII, set groupToPatch
        bool readASCII(IFstream&, const off_t);

        //- Read BINARY, set maxRegionId
        bool readBINARY(IFstream&, const off_t);

public:

    // Static Data

        //- The number of bytes in the STL binary header
        static const int headerSize = 80;

    // Static Member Functions

        //- Write "STL binary file" and number of triangles to stream
        static void writeHeaderBINARY(ostream&, unsigned int);

    // Constructors

        //- Read from file, filling in the information
        STLsurfaceFormatCore(const fileName&);

    // Destructor

        ~STLsurfaceFormatCore();

    // Member Functions

        //- File was detected to be STL binary
        bool binary() const
        {
            return binary_;
        }

        //- The number of facets is the number of regions
        label size() const
        {
            return regions_.size();
        }

        //- Flush all values
        void clear()
        {
            points_.clear();
            regions_.clear();
            groupToPatch_.clear();
        }

        //- Return full access to the faces
        pointField& points()
        {
            return points_;
        }

        //- Return full access to the regions
        List<label>& regions()
        {
            return regions_;
        }

        //- The max regionId when reading binary files
        label maxRegionId() const
        {
            return maxRegionId_;
        }

        //- The group/patch information when reading ASCII files
        const HashTable<label>& groupToPatch() const
        {
            return groupToPatch_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
