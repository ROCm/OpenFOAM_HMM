/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
     \\/     M anipulation  |
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
    base64Layer

Description
    An output filter layer to write base-64 encoded content.

    Base64 encoding accoding to RFC 4648 specification
    (https://tools.ietf.org/html/rfc4648#page-5).
    It is the obligation of the caller to avoid using normal output
    while the base-64 encoding layer is actively used.

SourceFiles
    base64Layer.C

\*---------------------------------------------------------------------------*/

#ifndef base64Layer_H
#define base64Layer_H

#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class base64Layer Declaration
\*---------------------------------------------------------------------------*/

class base64Layer
{
    // Private data

        //- The output stream for the layer
        std::ostream& os_;

        //- Buffer of characters to encode
        unsigned char group_[3];

        //- Current length of the encode buffer
        unsigned char groupLen_;

        //- Track if anything has been encoded.
        bool dirty_;


    // Private Member Functions

        inline unsigned char encode0() const;
        inline unsigned char encode1() const;
        inline unsigned char encode2() const;
        inline unsigned char encode3() const;

        //- Disallow default bitwise copy construct
        base64Layer(const base64Layer&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const base64Layer&) = delete;


protected:

    // Protected Member Functions

        //- Add a character to the group, outputting when the group is full.
        void add(char c);


public:

    // Constructors

        //- Construct and attach to an output stream
        base64Layer(std::ostream& os);


    //- Destructor
    ~base64Layer();


    // Member Functions

        //- The encoded length has 4 bytes out for every 3 bytes in.
        static std::size_t encodedLength(std::size_t n);


        //- Encode the character sequence, writing when possible.
        void write(const char* s, std::streamsize n);

        //- Restart a new encoding sequence.
        void reset();

        //- End the encoding sequence, padding the final characters with '='.
        //  Return false if no encoding was actually performed.
        bool close();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
