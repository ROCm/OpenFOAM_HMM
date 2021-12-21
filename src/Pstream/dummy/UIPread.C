/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

Description
    Read from UIPstream

\*---------------------------------------------------------------------------*/

#include "UIPstream.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::UIPstream::UIPstream
(
    const commsTypes commsType,
    const int fromProcNo,
    DynamicList<char>& receiveBuf,
    label& receiveBufPosition,
    const int tag,
    const label comm,
    const bool clearAtEnd,
    IOstreamOption::streamFormat fmt
)
:
    UPstream(commsType),
    Istream(fmt, IOstreamOption::currentVersion),
    fromProcNo_(fromProcNo),
    recvBuf_(receiveBuf),
    recvBufPos_(receiveBufPosition),
    tag_(tag),
    comm_(comm),
    clearAtEnd_(clearAtEnd),
    messageSize_(0)
{
    NotImplemented;
}


Foam::UIPstream::UIPstream(const int fromProcNo, PstreamBuffers& buffers)
:
    UPstream(buffers.commsType_),
    Istream(buffers.format_, IOstreamOption::currentVersion),
    fromProcNo_(fromProcNo),
    recvBuf_(buffers.recvBuf_[fromProcNo]),
    recvBufPos_(buffers.recvBufPos_[fromProcNo]),
    tag_(buffers.tag_),
    comm_(buffers.comm_),
    clearAtEnd_(true),
    messageSize_(0)
{
    NotImplemented;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::UIPstream::read
(
    const commsTypes commsType,
    const int fromProcNo,
    char* buf,
    const std::streamsize bufSize,
    const int tag,
    const label communicator
)
{
    NotImplemented;
    return 0;
}


// ************************************************************************* //
