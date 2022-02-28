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
    friend class UIPstreamBase;  // Access to recvBuf_, recvBufPos_


    // Private Data

        //- Track if sends are complete
        bool finishedSendsCalled_;

        //- Permit clear of individual receive buffer by external access
        bool allowClearRecv_;

        //- Buffer format (ascii | binary)
        const IOstreamOption::streamFormat format_;

        //- Communications type of this stream
        const UPstream::commsTypes commsType_;

        //- The transfer message type
        const int tag_;

        //- Communicator
        const label comm_;

        //- Send buffer. Size is nProcs()
        List<DynamicList<char>> sendBuf_;

        //- Receive buffer. Size is nProcs()
        List<DynamicList<char>> recvBuf_;

        //- Current read positions within recvBuf_. Size is nProcs()
        labelList recvBufPos_;


    // Private Member Functions

        //- Mark all sends as having been done.
        //  This will start receives (nonBlocking comms).
        void finalExchange(labelList& recvSizes, const bool wait);

        //- Mark sends as done.
        //  Only exchange sizes using the sendProcs/recvProcs subset
        //  (nonBlocking comms).
        void finalExchange
        (
            const labelUList& sendProcs,
            const labelUList& recvProcs,
            labelList& recvSizes,
            const bool wait
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


    // Sizing

        //- Number of ranks associated with PstreamBuffers
        label nProcs() const noexcept
        {
            return recvBufPos_.size();
        }

        //- Range of ranks indices associated with PstreamBuffers
        UPstream::rangeType allProcs() const noexcept
        {
            // Proc 0 -> nProcs (int value)
            return UPstream::rangeType(static_cast<int>(nProcs()));
        }


    // Queries

        //- True if finishedSends() or finishedNeighbourSends() has been called
        bool finished() const noexcept
        {
            return finishedSendsCalled_;
        }

        //- True if any (local) send buffers have data
        bool hasSendData() const;

        //- True if (local) send buffer has data on specified processor.
        bool hasSendData(const label proci) const;

        //- True if any (local) recv buffers have data.
        //- Must call finishedSends() or finishedNeighbourSends() first!
        bool hasRecvData() const;

        //- True if (local) recv buffer has data on specified processor.
        //- Must call finishedSends() or finishedNeighbourSends() first!
        bool hasRecvData(const label proci) const;

        //- Is clearStorage of individual receive buffer by external hooks
        //- allowed? (default: true)
        bool allowClearRecv() const noexcept
        {
            return allowClearRecv_;
        }


    // Edit

        //- Clear individual buffers and reset states.
        //  Does not clear individual buffer storage.
        void clear();

        //- Clear individual buffer storage and reset states.
        void clearStorage();

        //- Change allowClearRecv, return previous value
        bool allowClearRecv(bool on) noexcept;


    // Functions

        //- Mark sends as done
        //
        //  Non-blocking mode: populates receive buffers (all-to-all).
        //  \param wait wait for requests to complete (in nonBlocking mode)
        void finishedSends(const bool wait = true);

        //- Mark sends as done.
        //- Recovers the sizes (bytes) received.
        //
        //  Non-blocking mode: populates receive buffers (all-to-all).
        //  \param[out] recvSizes the sizes (bytes) received
        //  \param wait wait for requests to complete (in nonBlocking mode)
        //
        //  \warning currently only valid for nonBlocking comms.
        void finishedSends(labelList& recvSizes, const bool wait = true);

        //- Mark sends as done using subset of send/recv ranks
        //- to exchange data on.
        //
        //  Non-blocking mode: populates receive buffers.
        //  \param sendProcs ranks used for sends
        //  \param recvProcs ranks used for recvs
        //  \param wait wait for requests to complete (in nonBlocking mode)
        //
        //  \warning currently only valid for nonBlocking comms.
        void finishedSends
        (
            const labelUList& sendProcs,
            const labelUList& recvProcs,
            const bool wait = true
        );

        //- Mark sends as done using subset of send/recv ranks
        //- to exchange data on. Recovers the sizes (bytes) received.
        //
        //  Non-blocking mode: populates receive buffers.
        //
        //  \param sendProcs ranks used for sends
        //  \param recvProcs ranks used for recvs
        //  \param[out] recvSizes the sizes (bytes) received
        //  \param wait wait for requests to complete (in nonBlocking mode)
        //
        //  \warning currently only valid for nonBlocking comms.
        void finishedSends
        (
            const labelUList& sendProcs,
            const labelUList& recvProcs,
            labelList& recvSizes,
            const bool wait = true
        );

        //- Mark sends as done using subset of send/recv ranks
        //- and recover the sizes (bytes) received.
        //
        //  Non-blocking mode: populates receive buffers.
        //
        //  \param neighProcs ranks used for sends/recvs
        //  \param wait wait for requests to complete (in nonBlocking mode)
        //
        //  \warning currently only valid for nonBlocking comms.
        //  \note Same as finishedSends with identical sendProcs/recvProcs
        void finishedNeighbourSends
        (
            const labelUList& neighProcs,
            const bool wait = true
        )
        {
            finishedSends(neighProcs, neighProcs, wait);
        }

        //- Mark sends as done using subset of send/recv ranks
        //- and recover the sizes (bytes) received.
        //
        //  Non-blocking mode: it will populate receive buffers.
        //
        //  \param neighProcs ranks used for sends/recvs
        //  \param[out] recvSizes the sizes (bytes) received
        //  \param wait wait for requests to complete (in nonBlocking mode)
        //
        //  \warning currently only valid for nonBlocking mode.
        void finishedNeighbourSends
        (
            const labelUList& neighProcs,
            labelList& recvSizes,
            const bool wait = true
        )
        {
            finishedSends(neighProcs, neighProcs, recvSizes, wait);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
