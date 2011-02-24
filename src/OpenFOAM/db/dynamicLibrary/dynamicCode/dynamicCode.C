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

#include "dynamicCode.H"
#include "dynamicCodeContext.H"
#include "stringOps.H"
#include "IFstream.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "dictionary.H"
#include "dlLibraryTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::dynamicCode::allowSystemOperations
(
    Foam::debug::infoSwitch("allowSystemOperations", 0)
);


const Foam::word Foam::dynamicCode::codeTemplateEnvName
    = "FOAM_CODE_TEMPLATES";

const Foam::fileName Foam::dynamicCode::codeTemplateDirName
    = "codeTemplates/dynamicCode";


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::dynamicCode::checkSecurity
(
    const char* title,
    const dictionary& dict
)
{
    if (isAdministrator())
    {
        FatalIOErrorIn
        (
            title,
            dict
        )   << "This code should not be executed by someone with administrator"
            << " rights due to security reasons." << nl
            << "(it writes a shared library which then gets loaded "
            << "using dlopen)"
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::dynamicCode::copyAndFilter
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
            "dynamicCode::copyAndFilter()"
            " const"
        )   << "Failed opening for reading " << is.name()
            << exit(FatalError);
    }

    if (!os.good())
    {
        FatalErrorIn
        (
            "dynamicCode::copyAndFilter()"
            " const"
        )   << "Failed writing " << os.name()
            << exit(FatalError);
    }

    // Copy file while rewriting $VARS and ${VARS}
    string line;
    do
    {
        is.getLine(line);

        // expand according to mapping
        // expanding according to env variables might cause too many
        // surprises
        stringOps::inplaceExpand(line, mapping);

        os  << line.c_str() << nl;
    }
    while (is.good());
}


Foam::List<Foam::fileName>
Foam::dynamicCode::resolveTemplates(const UList<fileName>& names)
{
    // try to get template from FOAM_CODESTREAM_TEMPLATES
    const fileName templateDir(Foam::getEnv(codeTemplateEnvName));

    DynamicList<fileName> badFiles(names.size());
    List<fileName> resolved(names.size());

    label nResolved = 0;

    forAll(names, fileI)
    {
        const fileName& templateName = names[fileI];

        fileName file;
        if (!templateDir.empty() && isDir(templateDir))
        {
            file = templateDir/templateName;
            if (!isFile(file, false))
            {
                file.clear();
            }
        }

        // not found - fallback to ~OpenFOAM expansion
        if (file.empty())
        {
            file = findEtcFile(codeTemplateDirName/templateName);
        }

        if (file.empty())
        {
            badFiles.append(templateName);
        }
        else
        {
            resolved[nResolved++] = file;
        }
    }

    resolved.setSize(nResolved);

    if (!badFiles.empty())
    {
        FatalErrorIn
        (
            "dynamicCode::resolveTemplates(..)"
        )   << "Could not find the code template(s): "
            << badFiles << nl
            << "Under the $" << codeTemplateDirName
            << " directory or via via the ~OpenFOAM/"
            << codeTemplateDirName << " expansion"
            << exit(FatalError);
    }

    return resolved;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicCode::dynamicCode(const word& codeName)
:
    codeName_(codeName),
    codeDirName_(codeName)
{
    filterVars_.set("typeName", codeName_);
    filterVars_.set("SHA1sum", SHA1Digest().str());
}

Foam::dynamicCode::dynamicCode(const word& codeName, const word& codeDirName)
:
    codeName_(codeName),
    codeDirName_(codeDirName)
{
    filterVars_.set("typeName", codeName_);
    filterVars_.set("SHA1sum", SHA1Digest().str());
}


// Foam::dynamicCode::dynamicCode(const dynamicCode& dc)
// :
//     codeName_(dc.codeName_),
//     copyFiles_(dc.copyFiles_),
//     filesContents_(dc.filesContents_)
// {}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dynamicCode::clear()
{
    filterVars_.clear();
    filterFiles_.clear();
    createFiles_.clear();
    filterVars_.set("typeName", codeName_);
    filterVars_.set("SHA1sum", SHA1Digest().str());
}



void Foam::dynamicCode::addCreateFile
(
    const fileName& name,
    const string& contents
)
{
    createFiles_.append(fileAndContent(name, contents));
}


void Foam::dynamicCode::addFilterFile
(
    const fileName& name
)
{
    filterFiles_.append(name);
}


void Foam::dynamicCode::setFilterContext
(
    const dynamicCodeContext& context
)
{
    filterVars_.set("code", context.code());
    filterVars_.set("codeInclude", context.include());
    filterVars_.set("SHA1sum", context.sha1().str());
}


void Foam::dynamicCode::setFilterVariable
(
    const word& key,
    const string& value
)
{
    filterVars_.set(key, value);
}


Foam::fileName Foam::dynamicCode::codePath() const
{
    return stringOps::expand("$FOAM_CASE/dynamicCode/" + codeDirName_);
}


Foam::fileName Foam::dynamicCode::libPath() const
{
    return
    (
        stringOps::expand
        (
            "$FOAM_CASE/dynamicCode/platforms/$WM_OPTIONS/lib/lib"
        )
      + codeName_ + ".so"
    );
}


Foam::string Foam::dynamicCode::libTarget() const
{
    return "LIB = $(PWD)/../platforms/$(WM_OPTIONS)/lib/lib" + codeName_;
}



bool Foam::dynamicCode::copyFilesContents() const
{
    if (!allowSystemOperations)
    {
        FatalErrorIn
        (
            "dynamicCode::copyFilesContents(const fileName&) const"
        )   <<  "Loading a shared library using case-supplied code is not"
            << " enabled by default" << nl
            << "because of security issues. If you trust the code you can"
            << " enable this" << nl
            << "facility be adding to the InfoSwitches setting in the system"
            << " controlDict:" << nl << nl
            << "    allowSystemOperations 1" << nl << nl
            << "The system controlDict is either" << nl << nl
            << "    ~/.OpenFOAM/$WM_PROJECT_VERSION/controlDict" << nl << nl
            << "or" << nl << nl
            << "    $WM_PROJECT_DIR/etc/controlDict" << nl
            << endl
            << exit(FatalError);
    }

    List<fileName> resolvedFiles = resolveTemplates(filterFiles_);

    // Create dir
    const fileName outputDir = this->codePath();

    // Create dir
    mkDir(outputDir);

    // Copy/filter files
    forAll(resolvedFiles, fileI)
    {
        const fileName& srcFile = resolvedFiles[fileI];
        const fileName  dstFile(outputDir/srcFile.name());

        IFstream is(srcFile);
        //Info<< "Reading from " << is.name() << endl;
        if (!is.good())
        {
            FatalErrorIn
            (
                "dynamicCode::copyFilesContents(const fileName&)"
                " const"
            )   << "Failed opening " << srcFile
                << exit(FatalError);
        }

        OFstream os(dstFile);
        //Info<< "Writing to " << dstFile.name() << endl;
        if (!os.good())
        {
            FatalErrorIn
            (
                "dynamicCode::copyFilesContents(const fileName&)"
                " const"
            )   << "Failed writing " << dstFile
                << exit(FatalError);
        }

        // variables mapping
        copyAndFilter(is, os, filterVars_);
    }


    // Create files:
    forAll(createFiles_, fileI)
    {
        const fileName dstFile
        (
            outputDir/stringOps::expand(createFiles_[fileI].first())
        );

        mkDir(dstFile.path());
        OFstream os(dstFile);
        //Info<< "Writing to " << createFiles_[fileI].first() << endl;
        if (!os.good())
        {
            FatalErrorIn
            (
                "dynamicCode::copyFilesContents()"
                " const"
            )   << "Failed writing " << dstFile
                << exit(FatalError);
        }
        os  << createFiles_[fileI].second().c_str() << endl;
    }

    return true;
}


bool Foam::dynamicCode::wmakeLibso() const
{
    const Foam::string wmakeCmd("wmake libso " + this->codePath());
    Info<< "Invoking " << wmakeCmd << endl;

    if (Foam::system(wmakeCmd))
    {
        return false;
    }
    else
    {
        return true;
    }
}


bool Foam::dynamicCode::writeDigest
(
    const fileName& dirName,
    const SHA1Digest& sha1
) const
{
    mkDir(dirName);
    OFstream os(dirName/"SHA1Digest");
    os  << sha1;
    return os.good();
}


Foam::SHA1Digest Foam::dynamicCode::readDigest(const fileName& dirName) const
{
    IFstream is(dirName/"SHA1Digest");
    return SHA1Digest(is);
}


bool Foam::dynamicCode::upToDate(const SHA1Digest& sha1) const
{
    const fileName dirName = this->codePath();
    if (!exists(dirName/"SHA1Digest") || readDigest(dirName) != sha1)
    {
        writeDigest(dirName, sha1);
        return false;
    }
    else
    {
        return true;
    }
}


bool Foam::dynamicCode::upToDate(const dynamicCodeContext& context) const
{
    return upToDate(context.sha1());
}


// bool Foam::dynamicCode::openLibrary() const
// {
//     return dlLibraryTable::openLibrary(this->libPath(), false);
// }
//
//
// bool Foam::dynamicCode::closeLibrary() const
// {
//     return dlLibraryTable::closeLibrary(this->libPath(), false);
// }
//
//
// void* Foam::dynamicCode::findLibrary() const
// {
//     return dlLibraryTable::findLibrary(this->libPath());
// }



// bool Foam::dynamicCode::read(const dictionary& dict)
// {
//     dict.lookup("createFiles") >> createFiles_;
//     dict.lookup("filterFiles") >> filterFiles_;
//     dict.lookup("filterVariables") >> filterVariables_;
//
//     return true;
// }
//
//
// void Foam::dynamicCode::writeDict(Ostream& os) const
// {
//     os.writeKeyword("createFiles") << createFiles_
//         << token::END_STATEMENT << nl;
//
//     os.writeKeyword("filterFiles") << filterFiles_
//         << token::END_STATEMENT << nl;
//
//     os.writeKeyword("filterVariables") << filterVariables_
//         << token::END_STATEMENT << nl;
// }


// ************************************************************************* //
