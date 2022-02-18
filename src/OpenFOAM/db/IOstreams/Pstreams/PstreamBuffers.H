/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021-2022 OpenCFD Ltd.
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
    Foam::PstreamBuffers

Description
    Buffers for inter-processor communications streams (UOPstream, UIPstream).

    Use UOPstream to stream data into buffers, call finishedSends() to
    notify that data is in buffers and then use IUPstream to get data out
    of received buffers. Works with both blocking and nonBlocking. Does
    not make much sense with scheduled since there you would not need these
    explicit buffers.

    Example usage:
    \code
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        for (const int proci : Pstream::allProcs())
        {
            if (proci != Pstream::myProcNo())
            {
                someObject vals;

                UOPstream send(proci, pBufs);
                send << vals;
            }
        }

        pBufs.finishedSends();   // no-op for blocking

        for (const int proci : Pstream::allProcs())
        {
            if (proci != Pstream::myProcNo())
            {
                UIPstream recv(proci, pBufs);
                someObject vals(recv);
            }
        }
    \endcode

SourceFiles
    PstreamBuffers.C

\*---------------------------------------------------------------------------*/

#include "Pstream.H"

#ifndef Foam_PstreamBuffers_H
#define Foam_PstreamBuffers_H

#include "DynamicList.H"
#include "UPstream.H"
#include "IOstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class PstreamBuffers Declaration
\*---------------------------------------------------------------------------*/

class PstreamBuffers
{
    // Friendship
    friend class UOPstreamBase;  // Access to sendBuf_
    friend class UIPstreamBase;  // Access to recvBuf_ recvBufPos_;


    // Private Data

        //- Track if sends are complete
        bool finishedSendsCalled_;

        //- Buffer format (ascii | binary)
        const IOstreamOption::streamFormat format_;

        //- Communications type of this stream
        const UPstream::commsTypes commsType_;

        //- The transfer message type
        const int tag_;

        //- Communicator
        const label comm_;

        //- Send buffer
        List<DynamicList<char>> sendBuf_;

        //- Receive buffer
        List<DynamicList<char>> recvBuf_;

        //- Current read positions within recvBuf_
        labelList recvBufPos_;


    // Private Member Functions

        //- Mark all sends as having been done.
        //  This will start receives in non-blocking mode.
        void finalExchange(labelList& recvSizes, const bool block);

        //- Mark all sends as having been done.
        //  Only exchange sizes between sendProcs/recvProcs
        void finalExchange
        (
            const labelUList& sendProcs,
            const labelUList& recvProcs,
            labelList& recvSizes,
            const bool block
        );


public:

    // Constructors

        //- Construct given comms type, communication options, IO format
        explicit PstreamBuffers
        (
            const UPstream::commsTypes commsType,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        );


    //- Destructor - checks that all data have been consumed
    ~PstreamBuffers();


    // Member Functions

    // Access

        //- The associated buffer format (ascii | binary)
        IOstreamOption::streamFormat format() const noexcept
        {
            return format_;
        }

        //- The communications type of the stream
        UPstream::commsTypes commsType() const noexcept
        {
            return commsType_;
        }

        //- The transfer message type
        int tag() const noexcept
        {
            return tag_;
        }

        //- Communicator
        label comm() const noexcept
        {
            return comm_;
        }

        //- True if finishedSends has been called
        bool finished() const noexcept
        {
            return finishedSendsCalled_;
        }


    // Edit

        //- Reset (clear) individual buffers and reset state.
        //  Does not clear buffer storage
        void clear();


    // Functions

        //- Mark all sends as having been done.
        //  This will start receives in non-blocking mode.
        //  If block will wait for all transfers to finish
        //  (only relevant for nonBlocking mode)
        void finishedSends(const bool block = true);

        //- Mark all sends as having been done.
        //  Same as above but also returns sizes (bytes) received.
        //  \note currently only valid for non-blocking.
        void finishedSends(labelList& recvSizes, const bool block = true);

        //- Mark all sends as having been done, with limited
        //- send/recv processor connections.
        //  \note currently only valid for non-blocking.
        void finishedSends
        (
            const labelUList& sendProcs,
            const labelUList& recvProcs,
            const bool block = true
        );

        //- Mark all sends as having been done, with limited
        //- send/recv processor connections.
        //  Same as above but also returns sizes (bytes) received.
        //  \note currently only valid for non-blocking.
        void finishedSends
        (
            const labelUList& sendProcs,
            const labelUList& recvProcs,
            labelList& recvSizes,
            const bool block = true
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
