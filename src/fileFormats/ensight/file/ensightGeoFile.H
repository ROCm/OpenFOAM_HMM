/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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
    Foam::ensightGeoFile

Description
    Specialized Ensight output with extra geometry file header

\*---------------------------------------------------------------------------*/

#ifndef Foam_ensightGeoFile_H
#define Foam_ensightGeoFile_H

#include "ensightFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ensightGeoFile Declaration
\*---------------------------------------------------------------------------*/

class ensightGeoFile
:
    public ensightFile
{
    // Private Member Functions

        //- Initialize outputs the header information
        void init();

        //- No copy construct
        ensightGeoFile(const ensightGeoFile&) = delete;

        //- No copy assignment
        void operator=(const ensightGeoFile&) = delete;


public:

    // Static Functions

        //- Return a null ensightGeoFile
        static const ensightGeoFile& null()
        {
            return NullObjectRef<ensightGeoFile>();
        }


    // Constructors

        //- Construct from path-name.
        //  The path-name is adjusted for valid ensight file naming.
        explicit ensightGeoFile
        (
            const fileName& pathname,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        );

        //- Construct from path and name.
        //  Only the name portion is adjusted for valid ensight file naming.
        ensightGeoFile
        (
            const fileName& path,
            const fileName& name,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        );


    //- Destructor
    ~ensightGeoFile() = default;


    // Member Functions

    // Output

        //- Write keyword with trailing newline
        virtual Ostream& writeKeyword(const keyType& key);


    // Convenience Output Methods

        //- Begin a part (0-based index).
        using ensightFile::beginPart;

        //- Begin a "part" (0-based index), with a description.
        void beginPart(const label index, const string& description);

        //- Begin a "coordinates" block
        void beginCoordinates(const label npoints);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
