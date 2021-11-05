/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "IOobject.H"
#include "dictionary.H"
#include "foamVersion.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::IOstreamOption Foam::IOobject::parseHeader(const dictionary& headerDict)
{
    IOstreamOption streamOpt;  // == (ASCII, currentVersion)

    // Treat "version" as optional
    {
        token tok;
        if (headerDict.readIfPresent("version", tok))
        {
            streamOpt.version(tok);
        }
    }

    // Treat "format" as mandatory, could also as optional
    streamOpt.format(headerDict.get<word>("format"));

    headerClassName_ = headerDict.get<word>("class");

    const word headerObject(headerDict.get<word>("object"));

    // The "note" entry is optional
    headerDict.readIfPresent("note", note_);

    // The "arch" information may be missing
    string arch;
    if (headerDict.readIfPresent("arch", arch))
    {
        unsigned val = foamVersion::labelByteSize(arch);
        if (val) sizeofLabel_ = static_cast<unsigned char>(val);

        val = foamVersion::scalarByteSize(arch);
        if (val) sizeofScalar_ = static_cast<unsigned char>(val);
    }

    return streamOpt;
}


bool Foam::IOobject::readHeader(dictionary& headerDict, Istream& is)
{
    if (IOobject::debug)
    {
        InfoInFunction << "Reading header for file " << is.name() << endl;
    }

    // Check Istream not already bad
    if (!is.good())
    {
        if (rOpt_ == MUST_READ || rOpt_ == MUST_READ_IF_MODIFIED)
        {
            FatalIOErrorInFunction(is)
                << " stream not open for reading essential object from file "
                << is.relativeName()
                << exit(FatalIOError);
        }

        if (IOobject::debug)
        {
            SeriousIOErrorInFunction(is)
                << " stream not open for reading from file "
                << is.relativeName() << endl;
        }

        return false;
    }

    token firstToken(is);

    if (is.good() && firstToken.isWord("FoamFile"))
    {
        headerDict.read(is, false);  // Read sub-dictionary content

        IOstreamOption streamOpt = parseHeader(headerDict);

        is.format(streamOpt.format());
        is.version(streamOpt.version());
        is.setLabelByteSize(sizeofLabel_);
        is.setScalarByteSize(sizeofScalar_);
    }
    else
    {
        IOWarningInFunction(is)
            << "First token could not be read or is not 'FoamFile'"
            << nl << nl
            << "Check header is of the form:" << nl << endl;

        writeHeader(Info);

        return false;
    }

    // Check stream is still OK
    if (is.good())
    {
        objState_ = GOOD;
    }
    else
    {
        if (rOpt_ == MUST_READ || rOpt_ == MUST_READ_IF_MODIFIED)
        {
            FatalIOErrorInFunction(is)
                << " stream failure while reading header"
                << " on line " << is.lineNumber()
                << " of file " << is.relativeName()
                << " for essential object:" << name()
                << exit(FatalIOError);
        }

        if (IOobject::debug)
        {
            InfoInFunction
                << "Stream failure while reading header"
                << " on line " << is.lineNumber()
                << " of file " << is.relativeName() << endl;
        }

        objState_ = BAD;

        return false;
    }

    if (IOobject::debug)
    {
        Info<< " .... read" << endl;
    }

    return true;
}


bool Foam::IOobject::readHeader(Istream& is)
{
    dictionary headerDict;
    return IOobject::readHeader(headerDict, is);
}


// ************************************************************************* //
