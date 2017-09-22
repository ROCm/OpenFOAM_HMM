/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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
    Output to file stream, using a OSstream

SourceFiles
    OFstream.C

\*---------------------------------------------------------------------------*/

#ifndef OFstream_H
#define OFstream_H

#include "OSstream.H"
#include "fileName.H"
#include "className.H"

#include <fstream>
using std::ofstream;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class OFstreamAllocator Declaration
\*---------------------------------------------------------------------------*/

//- A std::ostream with the ability to handle compressed files
class OFstreamAllocator
{
protected:

    // Member Data

        //- The allocated stream pointer (ofstream or ogzstream).
        std::ostream* allocatedPtr_;


    // Constructors

        //- Construct from pathname
        OFstreamAllocator
        (
            const fileName& pathname,
            IOstream::compressionType compression=IOstream::UNCOMPRESSED,
            const bool append = false
        );


    //- Destructor
    ~OFstreamAllocator();


    // Protected Member Functions

        //- Delete the stream pointer
        void deallocate();

};


/*---------------------------------------------------------------------------*\
                          Class OFstream Declaration
\*---------------------------------------------------------------------------*/

class OFstream
:
    public OFstreamAllocator,
    public OSstream
{
public:

    // Declare name of the class and its debug switch
    ClassName("OFstream");


    // Constructors

        //- Construct from pathname
        OFstream
        (
            const fileName& pathname,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            compressionType compression=UNCOMPRESSED,
            const bool append = false
        );


    //- Destructor
    ~OFstream();


    // Member functions

      // Access

        //- Read/write access to the name of the stream
        using OSstream::name;


      // STL stream

        //- Access to underlying std::ostream
        virtual std::ostream& stdStream();

        //- Const access to underlying std::ostream
        virtual const std::ostream& stdStream() const;


      // Print

        //- Print description of IOstream to Ostream
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
