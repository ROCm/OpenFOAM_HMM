/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
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
    // Private Data

        int toProcNo_;

        DynamicList<char>& sendBuf_;

        const int tag_;

        const label comm_;

        const bool sendAtDestruct_;


    // Private Member Functions

        //- Prepare send buffer for count bytes of output,
        //- with specified alignment.
        inline void prepareBuffer(const size_t count, const size_t align);

        //- Generic write data to the send buffer, aligned by sizeof(T)
        template<class T>
        inline void writeToBuffer(const T& val);

        //- Write count bytes of data to the send buffer,
        //- using align byte alignment
        inline void writeToBuffer
        (
            const void* data,
            const size_t count,
            const size_t align
        );

        //- Add a single char to the send buffer. No alignment needed
        inline void putChar(const char c);

        //- Write string length and string content.
        //  The content includes the trailing nul char.
        inline void putString(const std::string& str);


public:

    // Constructors

        //- Construct given process index to write to using the given
        //- attached send buffer, optional communication characteristics
        //- and IO format
        UOPstream
        (
            const commsTypes commsType,
            const int toProcNo,
            DynamicList<char>& sendBuf,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm,
            const bool sendAtDestruct = true,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        );

        //- Construct given buffers
        UOPstream(const int toProcNo, PstreamBuffers& buffers);


    //- Destructor.
    ~UOPstream();


    // Member Functions

        // Inquiry

            //- Return flags of output stream
            virtual ios_base::fmtflags flags() const
            {
                return ios_base::fmtflags(0);
            }


        // Write Functions

            //- Write given buffer to given processor
            static bool write
            (
                const commsTypes commsType,
                const int toProcNo,
                const char* buf,
                const std::streamsize bufSize,
                const int tag = UPstream::msgType(),
                const label communicator = UPstream::worldComm
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
            virtual Ostream& write(const char* data, std::streamsize count);

            //- Low-level raw binary output.
            virtual Ostream& writeRaw(const char* data, std::streamsize count);

            //- Begin marker for low-level raw binary output.
            //  The count indicates the number of bytes for subsequent
            //  writeRaw calls.
            virtual bool beginRawWrite(std::streamsize count);

            //- End marker for low-level raw binary output.
            virtual bool endRawWrite()
            {
                return true;
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

            //- Get the current padding character
            //  \return previous padding character
            virtual char fill() const
            {
                return 0;
            }

            //- Set padding character for formatted field up to field width
            virtual char fill(const char)
            {
                return 0;
            }

            //- Get width of output field
            virtual int width() const
            {
                return 0;
            }

            //- Set width of output field
            //  \return previous width
            virtual int width(const int)
            {
                 return 0;
            }

            //- Get precision of output field
            virtual int precision() const
            {
                 return 0;
            }

            //- Set precision of output field
            //  \return old precision
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

        //- Print stream description to Ostream
        void print(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
