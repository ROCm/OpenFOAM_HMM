/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "PstreamBuffers.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{

    DynamicList<char> PstreamBuffers::nullBuf(0);
}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::PstreamBuffers::PstreamBuffers
(
    const UPstream::commsTypes commsType,
    const label tag,
    IOstream::streamFormat format,
    IOstream::versionNumber version
)
:
    commsType_(commsType),
    tag_(tag),
    format_(format),
    version_(version),
    sendBuf_(UPstream::nProcs()),
    recvBuf_(UPstream::nProcs()),
    finishedSendsCalled_(false)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::PstreamBuffers::finishedSends()
{
    finishedSendsCalled_ = true;

    if (commsType_ == UPstream::nonBlocking)
    {
        labelListList sizes;
        Pstream::exchange<DynamicList<char>, char>
        (
            sendBuf_,
            recvBuf_,
            sizes,
            tag_
        );
    }
}


void Foam::PstreamBuffers::finishedSends(labelListList& sizes)
{
    finishedSendsCalled_ = true;

    if (commsType_ == UPstream::nonBlocking)
    {
        labelListList sizes;
        labelListList send,recv;

        Pstream::exchange<DynamicList<char>, char>
        (
            sendBuf_,
            recvBuf_,
            sizes,
            tag_
        );
    }
    else
    {
        sizes.setSize(UPstream::nProcs());
        labelList& nsTransPs = sizes[UPstream::myProcNo()];
        nsTransPs.setSize(UPstream::nProcs());

        forAll(sendBuf_, procI)
        {
            nsTransPs[procI] = sendBuf_[procI].size();
        }

        // Send sizes across.
        label oldTag = UPstream::msgType();
        UPstream::msgType() = tag_;
        combineReduce(sizes, UPstream::listEq());
        UPstream::msgType() = oldTag;
    }
}


// ************************************************************************* //
