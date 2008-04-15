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
#include "IGeometryDescriptorImpl.H"
#include "LogEntry.H"
#include "Paths.H"

// Namespaces
#include "FoamXNameSpaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IGeometryDescriptorImpl::IGeometryDescriptorImpl
(
    const word& geometryName
)
:
    name_(geometryName),
    displayName_(geometryName),
    description_(geometryName)
{
    static const char* functionName =
        "FoamX::IGeometryDescriptorImpl::IGeometryDescriptorImpl"
        "(const char* geometryName)";

    LogEntry log(functionName, __FILE__, __LINE__);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IGeometryDescriptorImpl::IGeometryDescriptorImpl
(
    const word& geometryName,
    const dictionary& outerDict
)
:
    name_(geometryName)
{
    static const char* functionName =
        "FoamX::IGeometryDescriptorImpl::IGeometryDescriptorImpl"
        "(const char* geometryName, dictionary& outerDict)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Make sure the dictionary contains the geometry definition.
        if (!outerDict.found(geometryName))
        {
            throw FoamXError
            (
                E_FAIL,
                "Geometry definition dictionary '" + name_ + "' not found.",
                functionName,
                __FILE__, __LINE__
            );
        }

        const dictionary& geomDict(outerDict.subDict(geometryName));

        // Check that the dictionary has the minimum required information.
        if
        (
            !geomDict.found("displayName")
         || !geomDict.found("description"))
        {
            throw FoamXError
            (
                E_FAIL,
                "Malformed geometry definition dictionary for geometry type '"
               + name_ + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Get the patch field properties.
        geomDict.lookup("displayName")>> displayName_;
        geomDict.lookup("description")>> description_;
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IGeometryDescriptorImpl::~IGeometryDescriptorImpl()
{
    static const char* functionName =
        "FoamX::IGeometryDescriptorImpl::~IGeometryDescriptorImpl()";

    LogEntry log(functionName, __FILE__, __LINE__);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::IGeometryDescriptorImpl::name()
{
    static const char* functionName =
        "FoamX::IGeometryDescriptorImpl::name()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate string and return.
    return CORBA::string_dup(name_.c_str());
}

void FoamX::IGeometryDescriptorImpl::name(const char* newName)
{
    static const char* functionName =
        "FoamX::IGeometryDescriptorImpl::name(const char* newName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    name_ = newName;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::IGeometryDescriptorImpl::displayName()
{
    static const char* functionName =
        "FoamX::IGeometryDescriptorImpl::displayName()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate string and return.
    return CORBA::string_dup(displayName_.c_str());
}

void FoamX::IGeometryDescriptorImpl::displayName(const char* newName)
{
    static const char* functionName =
        "FoamX::IGeometryDescriptorImpl::displayName(const char* newName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    displayName_ = newName;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::IGeometryDescriptorImpl::description()
{
    static const char* functionName =
        "FoamX::IGeometryDescriptorImpl::description()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate string and return.
    return CORBA::string_dup(description_.c_str());
}

void FoamX::IGeometryDescriptorImpl::description(const char* newDescription)
{
    static const char* functionName =
        "FoamX::IGeometryDescriptorImpl::description"
        "(const char* newDescription)";

    LogEntry log(functionName, __FILE__, __LINE__);

    description_ = newDescription;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IGeometryDescriptorImpl::save(DictionaryWriter& dict)
{
    static const char* functionName =
        "FoamX::IGeometryDescriptorImpl::save(DictionaryWriter& dict)";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Save this patch information to the dictionary.
    dict.startSubDict(name_);

    dict.writeEntry("displayName", displayName_);
    dict.writeEntry("description", description_);

    dict.endSubDict();
}


// ************************************************************************* //
