/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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
    Foam::OFstream

Description
    Output to file stream, using an OSstream

SourceFiles
    OFstream.C

\*---------------------------------------------------------------------------*/

#ifndef OFstream_H
#define OFstream_H

#include "OSstream.H"
#include "className.H"
#include "fstreamPointer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class OFstream Declaration
\*---------------------------------------------------------------------------*/

class OFstream
:
    private Foam::ofstreamPointer,
    public OSstream
{
public:

    //- Declare type-name (with debug switch)
    ClassName("OFstream");


    // Constructors

        //- Construct a null output file stream.
        //  Behaves like \c /dev/null and is named accordingly
        explicit OFstream(std::nullptr_t);

        //- Construct from pathname
        explicit OFstream
        (
            const fileName& pathname,
            IOstreamOption streamOpt = IOstreamOption(),
            const bool append = false
        );

        //- Construct from pathname, format (version, compression)
        OFstream
        (
            const fileName& pathname,
            IOstreamOption::streamFormat fmt,
            IOstreamOption::versionNumber ver = IOstreamOption::currentVersion,
            IOstreamOption::compressionType cmp = IOstreamOption::UNCOMPRESSED,
            const bool append = false
        )
        :
            OFstream(pathname, IOstreamOption(fmt, ver, cmp), append)
        {}


    //- Destructor
    ~OFstream() = default;


    // Member Functions

        //- Read/write access to the name of the stream
        using OSstream::name;


    // STL stream

        //- Access to underlying std::ostream
        virtual std::ostream& stdStream();

        //- Const access to underlying std::ostream
        virtual const std::ostream& stdStream() const;

        //- Rewind the stream so that it may be written again.
        //- Reopens the file (truncation)
        virtual void rewind();


    // Print

        //- Print stream description
        void print(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Global predefined null output stream "/dev/null"
extern OFstream Snull;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
