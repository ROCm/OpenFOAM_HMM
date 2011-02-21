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
    ISstream& sourceStr,
    OSstream& destStr
) const
{
    if (!sourceStr.good())
    {
        FatalErrorIn
        (
            "codeStreamTools::copyAndExpand()"
            " const"
        )   << "Failed opening for reading " << sourceStr.name()
            << exit(FatalError);
    }

    if (!destStr.good())
    {
        FatalErrorIn
        (
            "codeStreamTools::copyAndExpand()"
            " const"
        )   << "Failed writing " << destStr.name() << exit(FatalError);
    }

    // Copy file whilst rewriting environment vars
    string line;
    do
    {
        sourceStr.getLine(line);
        line.expand(true, true);  // replace any envvars inside substitutions
        destStr<< line.c_str() << nl;
    }
    while (sourceStr.good());
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


Foam::codeStreamTools::codeStreamTools(const codeStreamTools& otf)
:
    name_(otf.name_),
    copyFiles_(otf.copyFiles_),
    filesContents_(otf.filesContents_)
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
            << " enabled by default" << endl
            << "because of security issues. If you trust the code you can"
            << " enable this" << endl
            << "facility be adding to the InfoSwitches setting in the system"
            << " controlDict:" << endl
            << endl
            << "    allowSystemOperations 1" << endl
            << endl
            << "The system controlDict is either" << endl
            << endl
            << "    ~/.OpenFOAM/$WM_PROJECT_VERSION/controlDict" << endl
            << endl
            << "or" << endl
            << endl
            << "    $WM_PROJECT_DIR/etc/controlDict" << endl
            << endl
            << exit(FatalError);
    }

    // Create dir
    mkDir(dir);

    //Info<< "Setting envvar typeName=" << name_ << endl;
    setEnv("typeName", name_, true);
    // Copy any template files
    forAll(copyFiles_, i)
    {
        const List<Pair<string> >& rules = copyFiles_[i].second();
        forAll(rules, j)
        {
            //Info<< "Setting envvar " << rules[j].first() << endl;
            setEnv(rules[j].first(), rules[j].second(), true);
        }

        const fileName sourceFile = fileName(copyFiles_[i].first()).expand();
        const fileName destFile = dir/sourceFile.name();

        IFstream sourceStr(sourceFile);
        //Info<< "Reading from " << sourceStr.name() << endl;
        if (!sourceStr.good())
        {
            FatalErrorIn
            (
                "codeStreamTools::copyFilesContents()"
                " const"
            )   << "Failed opening " << sourceFile << exit(FatalError);
        }

        OFstream destStr(destFile);
        //Info<< "Writing to " << destFile.name() << endl;
        if (!destStr.good())
        {
            FatalErrorIn
            (
                "codeStreamTools::copyFilesContents()"
                " const"
            )   << "Failed writing " << destFile << exit(FatalError);
        }

        copyAndExpand(sourceStr, destStr);
    }

    // Files that are always written:
    forAll(filesContents_, i)
    {
        fileName f = fileName(dir/filesContents_[i].first()).expand();

        mkDir(f.path());
        OFstream str(f);
        //Info<< "Writing to " << filesContents_[i].first() << endl;
        if (!str.good())
        {
            FatalErrorIn
            (
                "codeStreamTools::copyFilesContents()"
                " const"
            )   << "Failed writing " << f << exit(FatalError);
        }
        str << filesContents_[i].second().c_str() << endl;
    }
    return true;
}


Foam::string Foam::codeStreamTools::stripLeading(const string& s)
{
    label sz = s.size();
    if (sz > 0 && s[0] == '\n')
    {
        return s(1, sz-1);
    }
    else
    {
        return s;
    }
}


bool Foam::codeStreamTools::writeDigest
(
    const fileName& dir,
    const SHA1Digest& sha1
)
{
    OFstream str(dir/"SHA1Digest");
    str << sha1;
    return str.good();
}


Foam::SHA1Digest Foam::codeStreamTools::readDigest(const fileName& dir)
{
    IFstream str(dir/"SHA1Digest");
    return SHA1Digest(str);
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
