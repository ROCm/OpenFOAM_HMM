/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#ifndef ensightGeoFile_H
#define ensightGeoFile_H

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

        //- Initialize by outputting header information
        void initialize();

        //- No copy construct
        ensightGeoFile(const ensightGeoFile&) = delete;

        //- No copy assignment
        void operator=(const ensightGeoFile&) = delete;


public:

    // Static Functions

        //- Return a null ensightGeoFile
        inline static const ensightGeoFile& null()
        {
            return NullObjectRef<ensightGeoFile>();
        }


    // Constructors

        //- Construct from pathName.
        //  The entire pathName is checked for valid ensight naming.
        explicit ensightGeoFile
        (
            const fileName& pathname,
            IOstream::streamFormat format=IOstream::BINARY
        );

        //- Construct from path and name.
        //  Only the name portion is checked for valid ensight naming.
        ensightGeoFile
        (
            const fileName& path,
            const fileName& name,
            IOstream::streamFormat format=IOstream::BINARY
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
