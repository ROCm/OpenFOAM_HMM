/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#ifndef PstreamBuffers_H
#define PstreamBuffers_H

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
    friend class UOPstream;
    friend class UIPstream;

    // Private Data

        //- Communications type of this stream
        const UPstream::commsTypes commsType_;

        //- The transfer message type
        const int tag_;

        //- Communicator
        const label comm_;

        //- Buffer format (ascii | binary)
        const IOstreamOption::streamFormat format_;

        //- Send buffer
        List<DynamicList<char>> sendBuf_;

        //- Receive buffer
        List<DynamicList<char>> recvBuf_;

        //- Current read positions within recvBuf_
        labelList recvBufPos_;

        //- Track if sends are complete
        bool finishedSendsCalled_;


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

        //- Mark all sends as having been done.
        //  This will start receives in non-blocking mode.
        //  If block will wait for all transfers to finish
        //  (only relevant for nonBlocking mode)
        void finishedSends(const bool block = true);

        //- Mark all sends as having been done.
        //  Same as above but also returns sizes (bytes) received.
        //  \note currently only valid for non-blocking.
        void finishedSends(labelList& recvSizes, const bool block = true);

        //- Reset (clear) individual buffers and reset state.
        //  Does not clear buffer storage
        void clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
