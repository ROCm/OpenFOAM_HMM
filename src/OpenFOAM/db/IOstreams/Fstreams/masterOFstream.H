/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
    Foam::masterOFstream

Description
    Master-only drop-in replacement for OFstream.

SourceFiles
    masterOFstream.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_masterOFstream_H
#define Foam_masterOFstream_H

#include "StringStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class masterOFstream Declaration
\*---------------------------------------------------------------------------*/

class masterOFstream
:
    public OStringStream
{
    // Private Data

        //- The backend file name
        const fileName pathName_;

        //- Atomic file creation (ignored with append)
        const IOstreamOption::atomicType atomic_;

        //- Output file compression
        const IOstreamOption::compressionType compression_;

        //- Open file in append mode
        const IOstreamOption::appendType append_;

        //- Should file be written (on this processor)
        const bool valid_;


    // Private Member Functions

        //- Open file with checking and write append contents
        void checkWrite
        (
            const fileName& fName,
            const char* str,
            std::streamsize len
        );

        //- Open file with checking and write append contents
        void checkWrite(const fileName& fName, const std::string& s);

        //- Commit buffered information, including parallel gather as required
        void commit();


public:

    // Constructors

        //- Construct with specified atomic behaviour (with worldComm)
        //- from pathname, stream option, optional append
        masterOFstream
        (
            IOstreamOption::atomicType atomic,
            const fileName& pathname,
            IOstreamOption streamOpt = IOstreamOption(),
            IOstreamOption::appendType append = IOstreamOption::NON_APPEND,
            const bool valid = true
        );

        //- Construct (with worldComm)
        //- from pathname, stream option, optional append
        explicit masterOFstream
        (
            const fileName& pathname,
            IOstreamOption streamOpt = IOstreamOption(),
            IOstreamOption::appendType append = IOstreamOption::NON_APPEND,
            const bool valid = true
        )
        :
            masterOFstream
            (
                IOstreamOption::NON_ATOMIC,
                pathname,
                streamOpt,
                append,
                valid
            )
        {}


    //- Destructor - commits buffered information to file
    ~masterOFstream();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
