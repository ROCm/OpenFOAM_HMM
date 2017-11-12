/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::UOPstream

Description
    Output inter-processor communications stream operating on external
    buffer.

SourceFiles
    UOPstream.C

\*---------------------------------------------------------------------------*/

#include "Pstream.H"

#ifndef UOPstream_H
#define UOPstream_H

#include "UPstream.H"
#include "Ostream.H"
#include "DynamicList.H"
#include "PstreamBuffers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class UOPstream Declaration
\*---------------------------------------------------------------------------*/

class UOPstream
:
    public UPstream,
    public Ostream
{
    // Private data

        int toProcNo_;

        DynamicList<char>& sendBuf_;

        const int tag_;

        const label comm_;

        const bool sendAtDestruct_;


    // Private Member Functions

        //- Prepare buffer for count bytes of output at specified alignment.
        inline void prepareBuffer(const size_t count, const size_t align);

        //- Write data to the transfer buffer
        template<class T>
        inline void writeToBuffer(const T& val);

        //- Write a char to the transfer buffer
        inline void writeToBuffer(const char& c);

        //- Write count bytes of data to the transfer buffer
        //  using align byte alignment
        inline void writeToBuffer
        (
            const void* data,
            const size_t count,
            const size_t align
        );

        //- Write string length and content.
        //  The content includes the trailing nul char.
        inline void writeStringToBuffer(const std::string& str);


public:

    // Constructors

        //- Construct given process index to send to and optional buffer size,
        //  write format and IO version
        UOPstream
        (
            const commsTypes commsType,
            const int toProcNo,
            DynamicList<char>& sendBuf,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm,
            const bool sendAtDestruct = true,
            streamFormat format=BINARY,
            versionNumber version=currentVersion
        );

        //- Construct given buffers
        UOPstream(const int toProcNo, PstreamBuffers& buffers);


    //- Destructor
    ~UOPstream();


    // Member functions

        // Inquiry

            //- Return flags of output stream
            virtual ios_base::fmtflags flags() const
            {
                return ios_base::fmtflags(0);
            }


        // Write functions

            //- Write given buffer to given processor
            static bool write
            (
                const commsTypes commsType,
                const int toProcNo,
                const char* buf,
                const std::streamsize bufSize,
                const int tag = UPstream::msgType(),
                const label communicator = 0
            );

            //- Write token to stream or otherwise handle it.
            //  \return false if the token type was not handled by this method
            virtual bool write(const token& tok);

            //- Write single character. Whitespace is suppressed.
            virtual Ostream& write(const char c);

            //- Write the word-characters of a character string.
            //  Sends as a single char, or as word.
            virtual Ostream& write(const char* str);

            //- Write word
            virtual Ostream& write(const word& str);

            //- Write string
            virtual Ostream& write(const string& str);

            //- Write std::string surrounded by quotes.
            //  Optional write without quotes.
            virtual Ostream& writeQuoted
            (
                const std::string& str,
                const bool quoted=true
            );

            //- Write int32_t as a label
            virtual Ostream& write(const int32_t val);

            //- Write int64_t as a label
            virtual Ostream& write(const int64_t val);

            //- Write floatScalar
            virtual Ostream& write(const floatScalar val);

            //- Write doubleScalar
            virtual Ostream& write(const doubleScalar val);

            //- Write binary block with 8-byte alignment.
            virtual Ostream& write
            (
                const char* data,
                const std::streamsize count
            );

            //- Begin marker for low-level raw binary output.
            //  The count should indicate the number of bytes for subsequent
            //  writeRaw calls.
            virtual Ostream& beginRaw(const std::streamsize count);

            //- Low-level raw binary output.
            virtual Ostream& writeRaw
            (
                const char* data,
                const std::streamsize count
            );

            //- End marker for low-level raw binary output.
            virtual Ostream& endRaw()
            {
                return *this;
            }

            //- Add indentation characters
            virtual void indent()
            {}


        // Stream state functions

            //- Flush stream
            virtual void flush()
            {}

            //- Add newline and flush stream
            virtual void endl()
            {}

            //- Get width of output field
            virtual int width() const
            {
                return 0;
            }

            //- Set width of output field (and return old width)
            virtual int width(const int)
            {
                 return 0;
            }

            //- Get precision of output field
            virtual int precision() const
            {
                 return 0;
            }

            //- Set precision of output field (and return old precision)
            virtual int precision(const int)
            {
                 return 0;
            }


        // Edit

            //- Set flags of stream
            ios_base::fmtflags flags(const ios_base::fmtflags)
            {
                return ios_base::fmtflags(0);
            }


        // Print

            //- Print description of IOstream to Ostream
            void print(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
