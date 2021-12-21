/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::UIPstream

Description
    Input inter-processor communications stream operating on external
    buffer.

SourceFiles
    UIPstream.C

\*---------------------------------------------------------------------------*/

#include "Pstream.H"

#ifndef UIPstream_H
#define UIPstream_H

#include "UPstream.H"
#include "Istream.H"
#include "PstreamBuffers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class UIPstream Declaration
\*---------------------------------------------------------------------------*/

class UIPstream
:
    public UPstream,
    public Istream
{
    // Private Data

        int fromProcNo_;

        DynamicList<char>& recvBuf_;

        label& recvBufPos_;

        const int tag_;

        const label comm_;

        const bool clearAtEnd_;

        int messageSize_;


    // Private Member Functions

        //- Check buffer position against messageSize_ for EOF
        inline void checkEof();

        //- Prepare receive buffer by adjusting alignment
        inline void prepareBuffer(const size_t align);

        //- Read a T from the receive buffer
        template<class T>
        inline void readFromBuffer(T& val);

        //- Read count bytes of data from the receive buffer.
        //  Prior data alignment is done by prepareBuffer
        inline void readFromBuffer(void* data, const size_t count);

        //- Read string length and string content
        inline Istream& readString(std::string& str);


public:

    // Constructors

        //- Construct given process index to read from using the given
        //- attached receive buffer, optional communication characteristics
        //- and IO format
        UIPstream
        (
            const commsTypes commsType,
            const int fromProcNo,
            DynamicList<char>& receiveBuf,
            label& receiveBufPosition,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm,
            const bool clearAtEnd = false,   // destroy receiveBuf if at end
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        );

        //- Construct given buffers
        UIPstream(const int fromProcNo, PstreamBuffers& buffers);


    //- Destructor
    ~UIPstream();


    // Member Functions

        // Inquiry

            //- Return flags of output stream
            ios_base::fmtflags flags() const
            {
                return ios_base::fmtflags(0);
            }


        // Read Functions

            //- Read into given buffer from given processor
            //  \return the message size
            static label read
            (
                const commsTypes commsType,
                const int fromProcNo,
                char* buf,
                const std::streamsize bufSize,
                const int tag = UPstream::msgType(),
                const label communicator = UPstream::worldComm
            );

            //- Return next token from stream
            Istream& read(token& t);

            //- Read a character
            Istream& read(char& c);

            //- Read a word
            Istream& read(word& str);

            // Read a string
            Istream& read(string& str);

            //- Read a label
            Istream& read(label& val);

            //- Read a floatScalar
            Istream& read(floatScalar& val);

            //- Read a doubleScalar
            Istream& read(doubleScalar& val);

            //- Read binary block with 8-byte alignment.
            Istream& read(char* data, std::streamsize count);

            //- Low-level raw binary read
            Istream& readRaw(char* data, std::streamsize count);

            //- Start of low-level raw binary read
            bool beginRawRead();

            //- End of low-level raw binary read
            bool endRawRead()
            {
                return true;
            }

            //- Rewind the stream so that it may be read again
            void rewind();


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
