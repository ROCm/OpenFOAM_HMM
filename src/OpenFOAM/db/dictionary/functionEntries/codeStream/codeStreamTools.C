/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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

#include "codeStreamTools.H"
#include "stringOps.H"
#include "IFstream.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "dictionary.H"
#include "dlLibraryTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::codeStreamTools::allowSystemOperations
(
    Foam::debug::infoSwitch("allowSystemOperations", 0)
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::codeStreamTools::copyAndExpand
(
    ISstream& is,
    OSstream& os,
    const HashTable<string>& mapping
) const
{
    if (!is.good())
    {
        FatalErrorIn
        (
            "codeStreamTools::copyAndExpand()"
            " const"
        )   << "Failed opening for reading " << is.name()
            << exit(FatalError);
    }

    if (!os.good())
    {
        FatalErrorIn
        (
            "codeStreamTools::copyAndExpand()"
            " const"
        )   << "Failed writing " << os.name()
            << exit(FatalError);
    }

    // Copy file while rewriting $VARS and ${VARS}
    string line;
    do
    {
        is.getLine(line);

        // normal expansion according to mapping
        stringOps::inplaceExpand(line, mapping);

        // expand according to env variables
        stringOps::inplaceExpandEnv(line, true, true);

        os  << line.c_str() << nl;
    }
    while (is.good());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::codeStreamTools::codeStreamTools()
{}


Foam::codeStreamTools::codeStreamTools
(
    const word& name,
    const dictionary& dict
)
:
    name_(name)
{
    read(dict);
}


Foam::codeStreamTools::codeStreamTools
(
    const word& name,
    const List<fileAndVars>& copyFiles,
    const List<fileAndContent>& filesContents
)
:
    name_(name),
    copyFiles_(copyFiles),
    filesContents_(filesContents)
{}


Foam::codeStreamTools::codeStreamTools(const codeStreamTools& tools)
:
    name_(tools.name_),
    copyFiles_(tools.copyFiles_),
    filesContents_(tools.filesContents_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::codeStreamTools::copyFilesContents(const fileName& dir) const
{
    if (!allowSystemOperations)
    {
        FatalErrorIn
        (
            "codeStreamTools::copyFilesContents(const fileName&) const"
        )   <<  "Loading a shared library using case-supplied code is not"
            << " enabled by default" << nl
            << "because of security issues. If you trust the code you can"
            << " enable this" << nl
            << "facility be adding to the InfoSwitches setting in the system"
            << " controlDict" << nl
            << nl
            << "    allowSystemOperations 1" << nl << nl
            << "The system controlDict is either" << nl << nl
            << "    ~/.OpenFOAM/$WM_PROJECT_VERSION/controlDict" << nl << nl
            << "or" << nl << nl
            << "    $WM_PROJECT_DIR/etc/controlDict" << nl
            << endl
            << exit(FatalError);
    }

    // Create dir
    mkDir(dir);

    // Info<< "set mapping typeName=" << name_ << endl;
    // Copy any template files
    forAll(copyFiles_, i)
    {
        const fileName sourceFile(fileName(copyFiles_[i].file()).expand());
        const fileName destFile(dir/sourceFile.name());

        IFstream is(sourceFile);
        //Info<< "Reading from " << is.name() << endl;
        if (!is.good())
        {
            FatalErrorIn
            (
                "codeStreamTools::copyFilesContents(const fileName&)"
                " const"
            )   << "Failed opening " << sourceFile << exit(FatalError);
        }

        OFstream os(destFile);
        //Info<< "Writing to " << destFile.name() << endl;
        if (!os.good())
        {
            FatalErrorIn
            (
                "codeStreamTools::copyFilesContents(const fileName&)"
                " const"
            )   << "Failed writing " << destFile << exit(FatalError);
        }

        // variables mapping
        HashTable<string> mapping(copyFiles_[i]);
        mapping.set("typeName", name_);
        copyAndExpand(is, os, mapping);
    }


    // Files that are always written:
    forAll(filesContents_, i)
    {
        const fileName dstFile
        (
            fileName(dir/filesContents_[i].first()).expand()
        );

        mkDir(dstFile.path());
        OFstream os(dstFile);
        //Info<< "Writing to " << filesContents_[i].first() << endl;
        if (!os.good())
        {
            FatalErrorIn
            (
                "codeStreamTools::copyFilesContents()"
                " const"
            )   << "Failed writing " << dstFile << exit(FatalError);
        }
        os  << filesContents_[i].second().c_str() << endl;
    }

    return true;
}


bool Foam::codeStreamTools::writeDigest
(
    const fileName& dir,
    const SHA1Digest& sha1
)
{
    OFstream os(dir/"SHA1Digest");
    os  << sha1;
    return os.good();
}


Foam::SHA1Digest Foam::codeStreamTools::readDigest(const fileName& dir)
{
    IFstream is(dir/"SHA1Digest");
    return SHA1Digest(is);
}


bool Foam::codeStreamTools::upToDate
(
    const fileName& dir,
    const SHA1Digest& sha1
)
{
    if (!exists(dir/"SHA1Digest") || readDigest(dir) != sha1)
    {
        writeDigest(dir, sha1);
        return false;
    }
    else
    {
        return true;
    }
}


bool Foam::codeStreamTools::read(const dictionary& dict)
{
    dict.lookup("copyFiles") >> copyFiles_;
    dict.lookup("filesContents") >> filesContents_;

    return true;
}


void Foam::codeStreamTools::writeDict(Ostream& os) const
{
    os.writeKeyword("copyFiles") << copyFiles_ << token::END_STATEMENT << nl;
    os.writeKeyword("filesContents") << filesContents_ << token::END_STATEMENT
        << nl;
}


// ************************************************************************* //
