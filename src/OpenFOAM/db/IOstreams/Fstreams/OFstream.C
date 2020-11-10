/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "OFstream.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(OFstream, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::OFstream::OFstream
(
    std::nullptr_t
)
:
    Foam::ofstreamPointer(nullptr),
    OSstream(*(ofstreamPointer::get()), "/dev/null")
{
    setState(ofstreamPointer::get()->rdstate());
    setOpened();

    lineNumber_ = 1;
}


Foam::OFstream::OFstream
(
    const fileName& pathname,
    IOstreamOption streamOpt,
    const bool append
)
:
    Foam::ofstreamPointer(pathname, streamOpt.compression(), append),
    OSstream(*(ofstreamPointer::get()), pathname, streamOpt)
{
    setClosed();
    setState(ofstreamPointer::get()->rdstate());

    if (good())
    {
        setOpened();
    }
    else
    {
        setBad();
    }

    lineNumber_ = 1;

    if (debug)
    {
        if (pathname.empty())
        {
            InfoInFunction
                << "Cannot open empty file name"
                << Foam::endl;
        }

        if (!opened())
        {
            InfoInFunction
                << "Could not open file " << pathname
                << " for output\n" << info() << Foam::endl;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

std::ostream& Foam::OFstream::stdStream()
{
    std::ostream* ptr = ofstreamPointer::get();

    if (!ptr)
    {
        FatalErrorInFunction
            << "No stream allocated\n"
            << abort(FatalError);
    }

    return *ptr;
}


const std::ostream& Foam::OFstream::stdStream() const
{
    const std::ostream* ptr = ofstreamPointer::get();

    if (!ptr)
    {
        FatalErrorInFunction
            << "No stream allocated\n"
            << abort(FatalError);
    }

    return *ptr;
}


void Foam::OFstream::print(Ostream& os) const
{
    os  << "OFstream: ";
    OSstream::print(os);
}


// ************************************************************************* //
