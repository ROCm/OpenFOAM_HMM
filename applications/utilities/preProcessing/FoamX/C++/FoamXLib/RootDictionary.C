/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

// Foam header files
#include "IFstream.H"
#include "OSspecific.H"

// Project header files.
#include "FoamX.H"
#include "FoamXErrors.H"
#include "LogEntry.H"
#include "RootDictionary.H"

// Namespaces
#include "FoamXNameSpaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::RootDictionary::RootDictionary
(
    ITypeDescriptor_ptr typeDesc,
    const fileName& caseRoot,
    const fileName& caseName
)
:
    IDictionaryEntryImpl(typeDesc),
    caseRoot_(caseRoot),
    caseName_(caseName)
{
    static const char* functionName = 
        "FoamX::RootDictionary::RootDictionary"
        "(ITypeDescriptor_ptr, const fileName&, const fileName&)";

    LogEntry log(functionName, __FILE__, __LINE__);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::RootDictionary::~RootDictionary()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::fileName FoamX::RootDictionary::path
(
    const fileName& rootDir,
    const fileName& caseName,
    const fileName& dictPath
)
{
    // If path is relative prepend root and case
    if (dictPath[0] != '/')
    {
        return rootDir/caseName/dictPath;
    }
    else
    {
        return dictPath;
    }
}

Foam::fileName FoamX::RootDictionary::pathName()
{
    return 
        path(caseRoot_, caseName_, typeDescriptor_->dictionaryPath())
       /word(typeDescriptor_->name());
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::RootDictionary::load()
{
    fileName dictPathName = pathName();

    // See if the dictionary file exists.
    if (exists(dictPathName))
    {
        // Load the values from the file.
        load((IFstream(dictPathName)()));
    }
}

void FoamX::RootDictionary::load(Istream& is)
{
    static const char* functionName =
        "FoamX::RootDictionary::load(Istream&)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        log << "Loading file for "
            << CORBA::String_var(typeDescriptor_->path()) << endl;

        token firstToken(is);

        if
        (
            is.good()
         && firstToken.isWord()
         && firstToken.wordToken() == "FoamFile"
        )
        {
            // Read the FoamFile header
            dictionary headerDict(is);
        }
        else
        {
            is.putBack(firstToken);
        }

        IDictionaryEntryImpl::load(is);

        // Check state of Istream.
        is.check(functionName);
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::RootDictionary::save()
{
    static const char* functionName = "FoamX::RootDictionary::save()";

    // Overridden for root dictionary entries.
    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Make sure that this entry is a top-level compound type.
        fileName dictName(typeDescriptor_->name());

        if (!typeDescriptor_->isCompoundType() || !dictName.size())
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Invalid root dictionary object.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Get the dictionaries path and file name.
        fileName dictPath(typeDescriptor_->dictionaryPath());
        fileName dictPathName = dictPath/dictName;

        log << "Saving root dictionary " << dictName << "." << endl;

        if (dictPath[0] != '/')
        {
            DictionaryWriter dictWriter
            (
                caseRoot_,
                caseName_,
                dictPath,
                dictName
            );

            dictWriter.writeHeader
            (
                "FoamX Case Dictionary.",
                "dictionary"
            );

            writeEntries(dictWriter);
        }
        else
        {
            DictionaryWriter dictWriter(pathName());

            dictWriter.writeHeader
            (
                "Foam Dictionary.",
                "dictionary"
            );

            writeEntries(dictWriter);
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::RootDictionary::writeEntries(DictionaryWriter& dictWriter)
{
    static const char* functionName =
        "FoamX::RootDictionary::writeEntries(DictionaryWriter& dictWriter)";

    // Overridden for root dictionary entries.
    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Write the comment if required.
        const char* comment = typeDescriptor_->comment();
        if (strlen(comment) > 0)
        {
            dictWriter.writeComment(comment);
            dictWriter.writeEndl();
        }

        bool isDict = true;

        if (typeDescriptor_->type() != Type_Dictionary)
        {
            isDict = false;
        }

        // Save all sub-elements.
        for
        (
            DLList<IDictionaryEntryImpl*>::iterator iter =
                subElements_.begin();
            iter != subElements_.end();
            ++iter
        )
        {
            // Write each sub-element as a dictionary entry.
            iter()->save(dictWriter, isDict);
            dictWriter.writeEndl();
        }

        dictWriter.writeEndBar();
    }
    CATCH_ALL(functionName);
}


// ************************************************************************* //
