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

// FoamX header files.
#include "FoamX.H"
#include "FoamXErrors.H"
#include "IPatchDescriptorImpl.H"
#include "ITypeDescriptorImpl.H"
#include "LogEntry.H"
#include "Paths.H"

// Namespaces
#include "FoamXNameSpaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IPatchDescriptorImpl::IPatchDescriptorImpl(const word& patchName)
:
    name_(patchName),
    displayName_(patchName),
    description_(patchName)
{
    static const char* functionName =
        "FoamX::IPatchDescriptorImpl::IPatchDescriptorImpl"
        "(const char* patchName)";

    LogEntry log(functionName, __FILE__, __LINE__);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IPatchDescriptorImpl::IPatchDescriptorImpl
(
    const word& patchName,
    const dictionary& outerDict
)
:
    name_(patchName)
{
    static const char* functionName =
        "FoamX::IPatchDescriptorImpl::IPatchDescriptorImpl"
        "(const char* patchName, dictionary& outerDict)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Make sure the dictionary contains the patch field definition.
        if (!outerDict.found(patchName))
        {
            throw FoamXError
            (
                E_FAIL,
                "Patch definition dictionary '" + name_ + "' not found.",
                functionName,
                __FILE__, __LINE__
            );
        }

        const dictionary& patchDict(outerDict.subDict(patchName));

        // Check that the dictionary has the minimum required information.
        if
        (
            !patchDict.found("displayName")
         || !patchDict.found("description")
        )
        {
            throw FoamXError
            (
                E_FAIL,
                "Malformed patch definition dictionary for patch type '"
              + name_ + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Get the patch field properties.
        patchDict.lookup("displayName")>> displayName_;
        patchDict.lookup("description")>> description_;
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IPatchDescriptorImpl::~IPatchDescriptorImpl()
{
    static const char* functionName =
        "FoamX::IPatchDescriptorImpl::~IPatchDescriptorImpl()";

    LogEntry log(functionName, __FILE__, __LINE__);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::IPatchDescriptorImpl::name()
{
    static const char* functionName =
        "FoamX::IPatchDescriptorImpl::name()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate string and return.
    return CORBA::string_dup(name_.c_str());
}

void FoamX::IPatchDescriptorImpl::name(const char* newName)
{
    static const char* functionName =
        "FoamX::IPatchDescriptorImpl::name(const char* newName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    name_ = newName;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::IPatchDescriptorImpl::displayName()
{
    static const char* functionName =
        "FoamX::IPatchDescriptorImpl::displayName()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate string and return.
    return CORBA::string_dup(displayName_.c_str());
}

void FoamX::IPatchDescriptorImpl::displayName(const char* newName)
{
    static const char* functionName =
        "FoamX::IPatchDescriptorImpl::displayName(const char* newName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    displayName_ = newName;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::IPatchDescriptorImpl::description()
{
    static const char* functionName =
        "FoamX::IPatchDescriptorImpl::description()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate string and return.
    return CORBA::string_dup(description_.c_str());
}

void FoamX::IPatchDescriptorImpl::description(const char* newDescription)
{
    static const char* functionName =
        "FoamX::IPatchDescriptorImpl::description(const char* newDescription)";

    LogEntry log(functionName, __FILE__, __LINE__);

    description_ = newDescription;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPatchDescriptorImpl::save(DictionaryWriter& dict)
{
    static const char* functionName =
        "FoamX::IPatchDescriptorImpl::save(DictionaryWriter& dict)";

    LogEntry log(functionName, __FILE__, __LINE__);


    // Save this patch information to the dictionary.
    dict.startSubDict(name_);

    dict.writeEntry("displayName", displayName_);
    dict.writeEntry("description", description_);

    dict.endSubDict();
}


// ************************************************************************* //
