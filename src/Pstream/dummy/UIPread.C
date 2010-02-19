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

Description
    Read from UIPstream

\*---------------------------------------------------------------------------*/

#include "UIPstream.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::UIPstream::UIPstream
(
    const commsTypes commsType,
    const int fromProcNo,
    DynamicList<char>& externalBuf,
    label& externalBufPosition,
    const int tag,
    streamFormat format,
    versionNumber version
)
:
    UPstream(commsType),
    Istream(format, version),
    fromProcNo_(fromProcNo),
    externalBuf_(externalBuf),
    externalBufPosition_(externalBufPosition),
    tag_(tag),
    messageSize_(0)
{
    notImplemented
    (
        "UIPstream::UIPstream"
        "("
            "const commsTypes,"
            "const int fromProcNo,"
            "DynamicList<char>&,"
            "label&,"
            "const int tag,"
            "streamFormat, versionNumber"
        ")"
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::UIPstream::read
(
    const commsTypes commsType,
    const int fromProcNo,
    char* buf,
    const std::streamsize bufSize,
    const int tag
)
{
    notImplemented
    (
        "UIPstream::read"
        "("
            "const commsTypes,"
            "const int fromProcNo,"
            "char* buf,"
            "const label bufSize,"
            "const int tag"
        ")"
     );

     return 0;
}


// ************************************************************************* //
