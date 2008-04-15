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

// FoamX header files
#include "FoamX.H"
#include "FoamXErrors.H"
#include "IPatchPhysicalTypeDescriptorImpl.H"
#include "LogEntry.H"
#include "Paths.H"

// Namespaces
#include "FoamXNameSpaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IPatchPhysicalTypeDescriptorImpl::IPatchPhysicalTypeDescriptorImpl
(
    const word& patchPhysicalTypeName
)
:
    name_(patchPhysicalTypeName),
    displayName_(patchPhysicalTypeName),
    description_(patchPhysicalTypeName + " boundary condition"),
    patchType_("patch"),
    parentType_("")
{
    static const char* functionName =
        "FoamX::IPatchPhysicalTypeDescriptorImpl::"
        "IPatchPhysicalTypeDescriptorImpl"
        "(const char* patchPhysicalTypeName)";

    LogEntry log(functionName, __FILE__, __LINE__);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IPatchPhysicalTypeDescriptorImpl::IPatchPhysicalTypeDescriptorImpl
(
    const word& patchPhysicalTypeName,
    const wordList& fieldList
)
:
    name_(patchPhysicalTypeName),
    displayName_(patchPhysicalTypeName),
    description_(patchPhysicalTypeName + " boundary condition"),
    patchType_(patchPhysicalTypeName),
    parentType_("")
{
    static const char* functionName =
        "FoamX::IPatchPhysicalTypeDescriptorImpl::"
        "IPatchPhysicalTypeDescriptorImpl"
        "(const char* patchPhysicalTypeName, const wordList& fieldList)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Get list of fields
        patchFieldTypes_.length(fieldList.size());

        // Loop over all defined fields and determine the patch field types
        for (int nField = 0; nField < fieldList.size(); nField++)
        {
            patchFieldTypes_[nField].name = fieldList[nField].c_str();
            patchFieldTypes_[nField].value = patchPhysicalTypeName.c_str();
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool FoamX::IPatchPhysicalTypeDescriptorImpl::findParentPatchType
(
    word& patchType,
    const word& parentType,
    const dictionary& patchPhysicalTypesDict
) const
{
    if (patchPhysicalTypesDict.found(parentType))
    {
        const dictionary& parentPatchPhysicalTypeDict = 
            patchPhysicalTypesDict.subDict(parentType);

        if (parentPatchPhysicalTypeDict.found("patchType"))
        {
            parentPatchPhysicalTypeDict.lookup("patchType") >> patchType;
            return true;
        }
        else if (parentPatchPhysicalTypeDict.found("parentType"))
        {
            return findParentPatchType
            (
                patchType,
                parentPatchPhysicalTypeDict.lookup("parentType"),
                patchPhysicalTypesDict
            );
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


bool FoamX::IPatchPhysicalTypeDescriptorImpl::findPatchFieldType
(
    word& patchFieldType,
    const word& parentType,
    const dictionary& patchPhysicalTypesDict,
    const dictionary& patchFieldPhysicalTypes
) const
{
    if (patchFieldPhysicalTypes.found(parentType))
    {
        patchFieldPhysicalTypes.lookup(parentType) >> patchFieldType;
        return true;
    }
    else if
    (
        patchPhysicalTypesDict.found(parentType)
     && patchPhysicalTypesDict.subDict(parentType).found("parentType")
    )
    {
        return findPatchFieldType
        (
            patchFieldType,
            patchPhysicalTypesDict.subDict(parentType).lookup("parentType"),
            patchPhysicalTypesDict,
            patchFieldPhysicalTypes
        );
    }
    else
    {
        return false;
    }
}


FoamX::IPatchPhysicalTypeDescriptorImpl::IPatchPhysicalTypeDescriptorImpl
(
    const word& patchPhysicalTypeName,
    const dictionary& patchPhysicalTypeDict,
    const dictionary& patchPhysicalTypesDict,
    const dictionary& patchFieldsPhysicalTypesDict,
    const wordList& fieldList
)
:
    name_(patchPhysicalTypeName),
    displayName_(patchPhysicalTypeName),
    description_(patchPhysicalTypeName + " boundary condition"),
    patchType_("patch"),
    parentType_("")
{
    static const char* functionName =
        "FoamX::IPatchPhysicalTypeDescriptorImpl::"
        "IPatchPhysicalTypeDescriptorImpl"
        "(const char* patchPhysicalTypeName, "
        "dictionary& patchPhysicalTypeDict, "
        "dictionary& patchFieldsPhysicalTypesDict, "
        "const wordList& fieldList)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Get the boundary type properties
        if (patchPhysicalTypeDict.found("displayName"))
        {
            patchPhysicalTypeDict.lookup("displayName") >> displayName_;
        }

        if (patchPhysicalTypeDict.found("description"))
        {
            patchPhysicalTypeDict.lookup("description") >> description_;
        }

        if (patchPhysicalTypeDict.found("patchType"))
        {
            patchPhysicalTypeDict.lookup("patchType") >> patchType_;
        }

        // Get (optional) parent type
        if (patchPhysicalTypeDict.found("parentType"))
        {
            patchPhysicalTypeDict.lookup("parentType") >> parentType_;

            findParentPatchType
            (
                patchType_,
                parentType_,
                patchPhysicalTypesDict
            );
        }

        // Get list of fields
        patchFieldTypes_.length(fieldList.size());

        // Loop over all defined fields and determine the patch field types
        for (int nField=0; nField<fieldList.size(); nField++)
        {
            // Make sure there is an entry for this field
            // in patchFieldsPhysicalTypesDict
            if (!patchFieldsPhysicalTypesDict.found(fieldList[nField]))
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Field '" + fieldList[nField]
                  + "' does not exist in patchFieldsPhysicalTypes dictionary.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            const dictionary& patchFieldPhysicalTypes
                = patchFieldsPhysicalTypesDict.subDict(fieldList[nField]);

            word patchFieldType;

            // Find the patchPhysicalType for this field in 
            // patchFieldsPhysicalTypesDict hierachy
            if
            (
                findPatchFieldType
                (
                    patchFieldType,
                    patchPhysicalTypeName,
                    patchPhysicalTypesDict,
                    patchFieldPhysicalTypes
                )
            )
            {}
            else
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Boundary type '" + patchPhysicalTypeName
                  + "' does not exist in patchFieldsPhysicalTypes dictionary "
                  + "for field '" + fieldList[nField] + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            patchFieldTypes_[nField].name = fieldList[nField].c_str();
            patchFieldTypes_[nField].value = patchFieldType.c_str();
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IPatchPhysicalTypeDescriptorImpl::~IPatchPhysicalTypeDescriptorImpl()
{
    static const char* functionName =
        "FoamX::IPatchPhysicalTypeDescriptorImpl::"
        "~IPatchPhysicalTypeDescriptorImpl()";

    LogEntry log(functionName, __FILE__, __LINE__);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::IPatchPhysicalTypeDescriptorImpl::name()
{
    static const char* functionName =
        "FoamX::IPatchPhysicalTypeDescriptorImpl::name()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate string and return
    return CORBA::string_dup(name_.c_str());
}

void FoamX::IPatchPhysicalTypeDescriptorImpl::name(const char* newName)
{
    static const char* functionName =
        "FoamX::IPatchPhysicalTypeDescriptorImpl::name(const char* newName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    name_ = newName;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::IPatchPhysicalTypeDescriptorImpl::displayName()
{
    static const char* functionName =
        "FoamX::IPatchPhysicalTypeDescriptorImpl::displayName()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate string and return
    return CORBA::string_dup(displayName_.c_str());
}

void FoamX::IPatchPhysicalTypeDescriptorImpl::displayName(const char* newName)
{
    static const char* functionName =
        "FoamX::IPatchPhysicalTypeDescriptorImpl::"
        "displayName(const char* newName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    displayName_ = newName;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::IPatchPhysicalTypeDescriptorImpl::description()
{
    static const char* functionName =
        "FoamX::IPatchPhysicalTypeDescriptorImpl::description()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate string and return
    return CORBA::string_dup(description_.c_str());
}

void FoamX::IPatchPhysicalTypeDescriptorImpl::description
(
    const char* newDescription
)
{
    static const char* functionName =
        "FoamX::IPatchPhysicalTypeDescriptorImpl::description"
        "(const char* newDescription)";

    LogEntry log(functionName, __FILE__, __LINE__);

    description_ = newDescription;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::IPatchPhysicalTypeDescriptorImpl::patchType()
{
    static const char* functionName =
        "FoamX::IPatchPhysicalTypeDescriptorImpl::patchType()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate string and return
    return CORBA::string_dup(patchType_.c_str());
}

void FoamX::IPatchPhysicalTypeDescriptorImpl::patchType(const char* newName)
{
    static const char* functionName =
        "FoamX::IPatchPhysicalTypeDescriptorImpl::patchType"
        "(const char* newName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    patchType_ = newName;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::IPatchPhysicalTypeDescriptorImpl::parentType()
{
    static const char* functionName =
        "FoamX::IPatchPhysicalTypeDescriptorImpl::parentType()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate string and return
    return CORBA::string_dup(parentType_.c_str());
}

void FoamX::IPatchPhysicalTypeDescriptorImpl::parentType(const char* newName)
{
    static const char* functionName =
        "FoamX::IPatchPhysicalTypeDescriptorImpl::parentType"
        "(const char* newName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    parentType_ = newName;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::StringPairList*
FoamX::IPatchPhysicalTypeDescriptorImpl::patchFieldTypes()
{
    static const char* functionName =
        "FoamX::IPatchPhysicalTypeDescriptorImpl::patchFieldTypes()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate list and return
    return new StringPairList(patchFieldTypes_);
}

void FoamX::IPatchPhysicalTypeDescriptorImpl::patchFieldTypes
(
    const FoamXServer::StringPairList& newList
)
{
    static const char* functionName =
        "FoamX::IPatchPhysicalTypeDescriptorImpl::patchFieldTypes";

    LogEntry log(functionName, __FILE__, __LINE__);

    patchFieldTypes_ = newList;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IPatchPhysicalTypeDescriptorImpl::save
(
    DictionaryWriter& dictWriter
)
{
    static const char* functionName =
        "FoamX::IPatchPhysicalTypeDescriptorImpl::save"
        "(DictionaryWriter& dictWriter)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {

        // Save the boundary type information to a sub-dictionary
        dictWriter.startSubDict(name_);

        dictWriter.writeEntry("displayName", displayName_);
        dictWriter.writeEntry("description", description_);
        dictWriter.writeEntry("patchType", patchType_);

        // Only write the parent type if it has been set
        if (parentType_.size() > 0)
        {
            dictWriter.writeEntry("parentType", parentType_);
        }

        // Write the patch field types
        for (unsigned int i=0; i<patchFieldTypes_.length(); i++)
        {
            dictWriter.writeEntry
            (
                word(patchFieldTypes_[i].name),
                word(patchFieldTypes_[i].value)
            );
        }

        dictWriter.endSubDict();
    }
    CATCH_ALL(functionName);
}


// ************************************************************************* //
