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
#include "polyPatch.H"
#include "OSspecific.H"

// FoamX header files
#include "FoamX.H"
#include "FoamXErrors.H"
#include "ITypeDescriptorImpl.H"
#include "IApplicationImpl.H"
#include "IPropertiesImpl.H"
#include "IGeometricFieldDescriptorImpl.H"
#include "IPatchPhysicalTypeDescriptorImpl.H"
#include "DictionaryWriter.H"
#include "LogEntry.H"
#include "Paths.H"

// Namespaces
#include "FoamXNameSpaces.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IApplicationImpl::IApplicationImpl
(
    const FoamXServer::ApplicationDescriptor& appDesc,
    const IPropertiesImpl& foamSystemProperties
)
:
    name_(appDesc.name),
    description_(appDesc.name),
    category_(appDesc.category),
    appCfgPath_(fileName(appDesc.path)/"FoamX"),
    systemClass_(appDesc.systemClass),
    foamSystemProperties_(foamSystemProperties)
{
    static const char* functionName =
        "FoamX::IApplicationImpl::IApplicationImpl"
        "(const char* className, bool systemClass)";

    LogEntry log(functionName, __FILE__, __LINE__);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IApplicationImpl::~IApplicationImpl()
{
    static const char* functionName =
        "FoamX::IApplicationImpl::~IApplicationImpl()";

    LogEntry log(functionName, __FILE__, __LINE__);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IApplicationImpl::validate()
{
    static const char* functionName =
        "FoamX::IApplicationImpl::validate()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IApplicationImpl::load()
{
    static const char* functionName =
        "FoamX::IApplicationImpl::load()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Make sure that the definition dictionary exists
        if (!exists(appCfgPathName()))
        {
            throw FoamXError
            (
                E_FAIL,
                "Application Class definition file '"
               + appCfgPathName() + "' could not be found.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Set the FOAMX_APP_CFG_PATH environment variable
        // to the application configuration path
        Foam::setEnv("FOAMX_APP_CFG_PATH", appCfgPath_, true);

        // Open the application class definition dictionary
        dictionary classDict((IFstream(appCfgPathName())()));

        // Make sure we have all of the required entries
        if
        (
            !classDict.found("description")
         || !classDict.found("dictionaries")
         || !classDict.found("fields")
         || !classDict.found("patchPhysicalTypes")
         || !classDict.found("patchFieldsPhysicalTypes")
        )
        {
            throw FoamXError
            (
                E_FAIL,
                "Invalid application class configuration dictionary '"
              + appCfgPathName() + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Get application class info
        classDict.lookup("description") >> description_;

        if (classDict.found("modules"))
        {
            moduleList_.read(classDict.lookup("modules"));
        }
        else
        {
            moduleList_.append("FoamX.Modules.CaseEditor.CaseEditorModule");
        }

        const dictionary& dictionariesDict = classDict.subDict("dictionaries");

        // Initialise the dictionary objects
        for
        (
            dictionary::const_iterator iter = dictionariesDict.begin();
            iter != dictionariesDict.end();
            ++iter
        )
        {
            const word& dictionaryName = iter().keyword();

            // Set/correct the application name
            if
            (
                dictionaryName == "controlDict"
             && iter().isDict()
             && iter().dict().found("entries")
             && iter().dict().subDict("entries").found("application")
            )
            {
                dictionary& applicationNameDict = const_cast<dictionary&>
                (
                    iter().dict().subDict("entries").subDict("application")
                );

                applicationNameDict.remove("default");
                applicationNameDict.add("default", name_);
            }

            // See if we already have a dictionary with this name
            if (dictTypeDescriptorMap_.found(dictionaryName))
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Duplicate dictionary name '" + dictionaryName + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Create and initialise a new TypeDescriptor object for this
            // dictionary
            ITypeDescriptorImpl* pTypeDescriptor = new ITypeDescriptorImpl
            (
                dictionaryName,
                appCfgPath_,
                iter(),
                foamSystemProperties_.foamTypesDict()
            );

            if (pTypeDescriptor == NULL)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Failed to create Dictionary TypeDescriptor object for "
                    "dictionary '" + dictionaryName + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            dictTypeDescriptorMap_.append(dictionaryName, pTypeDescriptor);
        }

        /* This does work but should only be used if it is considered useful
        if (classDict.found("types"))
        {
            (dictionary&)foamSystemProperties_.foamTypesDict()
                += classDict.subDict("types");
        }
        */

        const dictionary& fieldsDict = classDict.subDict("fields");

        // Initialise the field objects
        for
        (
            dictionary::const_iterator iter = fieldsDict.begin();
            iter != fieldsDict.end();
            ++iter
        )
        {
            const word& fieldName = iter().keyword();

            // See if we already have a field with this name
            if (fieldDescriptorMap_.found(fieldName))
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Duplicate field name '" + fieldName + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Create and initialise a new FieldDescriptor object for field
            IGeometricFieldDescriptorImpl* pFieldDescriptor = new
            IGeometricFieldDescriptorImpl
            (
                fieldName,
                fieldsDict.subDict(fieldName),
                foamSystemProperties_.foamTypes(),
                foamSystemProperties_.geometryDescriptors()
            );

            if (pFieldDescriptor == NULL)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Failed to create FieldDescriptor object for field '"
                   + fieldName + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            fieldDescriptorMap_.append(fieldName, pFieldDescriptor);
        }

        wordList fieldList = fieldDescriptorMap_.toc();


        wordList constranedPatchTypes = polyPatch::constraintTypes();

        forAll(constranedPatchTypes, i)
        {
            // See if we already have a boundary type with this name
            if (patchPhysicalTypeDescriptorMap_.found(constranedPatchTypes[i]))
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Duplicate boundary type name '"
                  + constranedPatchTypes[i] + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Create and initialise a new PatchPhysicalTypeDescriptor object
            // for this field
            IPatchPhysicalTypeDescriptorImpl* pBoundaryDescriptor = new
            IPatchPhysicalTypeDescriptorImpl
            (
                constranedPatchTypes[i],
                fieldList
            );

            if (pBoundaryDescriptor == NULL)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Failed to create PatchPhysicalTypeDescriptor object for "
                    "boundary type '" + constranedPatchTypes[i] + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            patchPhysicalTypeDescriptorMap_.append
            (
                constranedPatchTypes[i],
                pBoundaryDescriptor
            );
        }

        if (classDict.found("patchFieldTypes"))
        {
            const_cast<IPropertiesImpl&>(foamSystemProperties_).addPatchFields
            (
                classDict.subDict("patchFieldTypes")
            );
        }

        const dictionary& patchPhysicalTypesDict
        (
            classDict.subDict("patchPhysicalTypes")
        );

        const dictionary& patchFieldsPhysicalTypesDict
        (
            classDict.subDict("patchFieldsPhysicalTypes")
        );

        // Initialise the boundary type information
        for
        (
            dictionary::const_iterator iter = patchPhysicalTypesDict.begin();
            iter != patchPhysicalTypesDict.end();
            ++iter
        )
        {
            const word& patchPhysicalTypeName = iter().keyword();

            // See if we already have a boundary type with this name
            if (patchPhysicalTypeDescriptorMap_.found(patchPhysicalTypeName))
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Duplicate boundary type name '"
                  + patchPhysicalTypeName + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Create and initialise a new PatchPhysicalTypeDescriptor object
            // for this field
            IPatchPhysicalTypeDescriptorImpl* pBoundaryDescriptor = new
            IPatchPhysicalTypeDescriptorImpl
            (
                patchPhysicalTypeName,
                iter().dict(),
                patchPhysicalTypesDict,
                patchFieldsPhysicalTypesDict,
                fieldList
            );

            if (pBoundaryDescriptor == NULL)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Failed to create PatchPhysicalTypeDescriptor object for "
                    "boundary type '" + patchPhysicalTypeName + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            patchPhysicalTypeDescriptorMap_.append
            (
                patchPhysicalTypeName,
                pBoundaryDescriptor
            );
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IApplicationImpl::save()
{
    static const char* functionName =
        "FoamX::IApplicationImpl::save()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Cannot save system application classes
        if (systemClass_)
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Cannot save system application class definition file '"
               + appCfgPathName() + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Write the application class data
        DictionaryWriter dict(appCfgPathName());

        dict.writeHeader
        (
            "FoamX Application Class Configuration File for " + name_,
            "dictionary"
        );

        // Write execution info
        dict.writeSectionHeader("Execution");
        dict.writeEntry("description", description_);
        dict.writeEndl();
        dict.writeEntry("arguments", arguments_);

        // Write module names
        dict.writeSectionHeader("FoamX Modules");
        dict.writeEntry("modules", moduleList_);

        // Write dictionary names
        dict.writeSectionHeader("Dictionaries");
        wordList dictNames;
        dictNames.setSize(dictTypeDescriptorMap_.size());
        label i = 0;
        for
        (
            Dictionary<ITypeDescriptorImpl>::iterator iter =
                dictTypeDescriptorMap_.begin();
            iter != dictTypeDescriptorMap_.end();
            ++iter
        )
        {
            CORBA::String_var name = iter().name();
            dictNames[i++] = word(name);
        }

        // Use the top-level dictionary names instead of the keys
        dict.writeEntry("dictionaries", dictNames);

        // Write field definitions
        dict.writeSectionHeader("Fields");
        dict.startSubDict("fields");
        i = 0;

        for
        (
            Dictionary<IGeometricFieldDescriptorImpl>::iterator iter
                = fieldDescriptorMap_.begin();
            iter != fieldDescriptorMap_.end();
            ++iter
        )
        {
            iter().save(dict);

            if (i++ < fieldDescriptorMap_.size()-1)
            {
                dict.writeEndl();
            }
        }
        dict.endSubDict();

        // Write boundary types
        dict.writeSectionHeader("Boundary Types");
        dict.startSubDict("patchPhysicalTypes");
        i = 0;

        for
        (
            Dictionary<IPatchPhysicalTypeDescriptorImpl>::iterator iter
                = patchPhysicalTypeDescriptorMap_.begin();
            iter != patchPhysicalTypeDescriptorMap_.end();
            ++iter
        )
        {
            iter().save(dict);

            if (i++ < patchPhysicalTypeDescriptorMap_.size()-1)
            {
                dict.writeEndl();
            }
        }
        dict.endSubDict();

        dict.writeEndl();
        dict.writeEndBar();


        // Save the dictionary type definitions
        for
        (
            Dictionary<ITypeDescriptorImpl>::iterator iter =
                dictTypeDescriptorMap_.begin() ;
            iter != dictTypeDescriptorMap_.end() ;
            ++iter
        )
        {
            const word& dictKey = iter().name();

            ITypeDescriptorImpl* pDictType = dictTypeDescriptorMap_[dictKey];

            // Create a dictionary writer object for this dictionary

            // Use the dictionary name, not the key
            CORBA::String_var dictName = pDictType->name();

            fileName dictFilePath = appCfgPath_/word(dictName) + ".cfg";

            // Save dictionary type description
            DictionaryWriter dictWriter(dictFilePath);

            dictWriter.writeHeader
            (
                "FoamX Application Class Configuration File.",
                "dictionary"
            );

            // Top-level type descriptor
            pDictType->save(dictWriter, true);

            dictWriter.writeEndl();
            dictWriter.writeEndBar();
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::IApplicationImpl::name()
{
    static const char* functionName =
        "FoamX::IApplicationImpl::name()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate string and return
    return CORBA::string_dup(name_.c_str());
}

void FoamX::IApplicationImpl::name(const char* newName)
{
    static const char* functionName =
        "FoamX::IApplicationImpl::name(const char* newName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    name_ = newName;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::IApplicationImpl::description()
{
    static const char* functionName =
        "FoamX::IApplicationImpl::description()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate string and return
    return CORBA::string_dup(description_.c_str());
}

void FoamX::IApplicationImpl::description(const char* newDescription)
{
    static const char* functionName =
        "FoamX::IApplicationImpl::description(const char* newDescription)";

    LogEntry log(functionName, __FILE__, __LINE__);

    description_ = newDescription;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::IApplicationImpl::category()
{
    static const char* functionName =
        "FoamX::IApplicationImpl::category()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate string and return
    return CORBA::string_dup(category_.c_str());
}

void FoamX::IApplicationImpl::category(const char* newCategory)
{
    static const char* functionName =
        "FoamX::IApplicationImpl::category(const char* newCategory)";

    LogEntry log(functionName, __FILE__, __LINE__);

    category_ = newCategory;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::StringList* FoamX::IApplicationImpl::modules()
{
    static const char* functionName =
        "FoamX::IApplicationImpl::modules()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate list and return
    return new FoamXServer::StringList(moduleList_);
}

void FoamX::IApplicationImpl::modules
(
    const FoamXServer::StringList& newList
)
{
    static const char* functionName =
        "FoamX::IApplicationImpl::modules"
        "(const FoamXServer::StringList& newList)";

    LogEntry log(functionName, __FILE__, __LINE__);

    moduleList_ = newList;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

CORBA::Boolean FoamX::IApplicationImpl::systemClass()
{
    static const char* functionName =
        "FoamX::IApplicationImpl::systemClass()";

    LogEntry log(functionName, __FILE__, __LINE__);

    return systemClass_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::StringList* FoamX::IApplicationImpl::fields()
{
    static const char* functionName =
        "FoamX::IApplicationImpl::fields()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate list and return
    return new FoamXServer::StringList
    (
        FoamXWordList(static_cast<const wordList&>(fieldDescriptorMap_.toc()))
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IApplicationImpl::getField
(
    const char* fieldName,
    IGeometricFieldDescriptor_out fieldDescriptor
)
{
    static const char* functionName =
        "FoamX::IApplicationImpl::getField"
        "(const char* fieldName, IGeometricFieldDescriptor_out fieldDescriptor)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {

        // See if we have this field's type descriptor cached
        if (!fieldDescriptorMap_.found(fieldName))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid field name.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Return a reference to the FieldDescriptor object
        fieldDescriptor = fieldDescriptorMap_[fieldName]->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void  FoamX::IApplicationImpl::findField
(
    const char* fieldName,
    IGeometricFieldDescriptor_out fieldDescriptor
)
{
    static const char* functionName =
        "FoamX::IApplicationImpl::findField"
        "(const char* fieldName, IGeometricFieldDescriptor_out fieldDescriptor)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {

        // Loop over all known field descriptors and find the one with the
        // specified name
        for
        (
            ObjRefDictionary<IGeometricFieldDescriptorImpl>::iterator iter =
                fieldDescriptorMap_.begin() ;
            iter != fieldDescriptorMap_.end() ;
            ++iter
        )
        {
            // Check for matching name
            CORBA::String_var fieldDescName = iter().name();
            if (strcmp(fieldDescName, fieldName) == 0)
            {
                fieldDescriptor = iter()._this();
                break;
            }
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IApplicationImpl::addField
(
    const char* fieldName,
    IGeometricFieldDescriptor_out fieldDescriptor
)
{
    static const char* functionName =
        "FoamX::IApplicationImpl::addField(const char* fieldName, "
        "IGeometricFieldDescriptor_out fieldDescriptor)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // See if we already have a field with this name
        if (fieldDescriptorMap_.found(fieldName))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid field name.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Create and initialise a new FieldDescriptor
        IGeometricFieldDescriptorImpl* pFieldDescriptor =
            new IGeometricFieldDescriptorImpl(fieldName);

        if (pFieldDescriptor == NULL)
        {
            throw FoamXError
            (
                E_FAIL,
                "Couldn't create FieldDescriptor object.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Add to field map
        fieldDescriptorMap_.append(fieldName, pFieldDescriptor);

        // Return a reference to the FieldDescriptor object
        fieldDescriptor = pFieldDescriptor->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IApplicationImpl::deleteField(const char* fieldName)
{
    static const char* functionName =
        "FoamX::IApplicationImpl::deleteField(const char* fieldName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {

        // See if we have a field with this name
        if (!fieldDescriptorMap_.found(fieldName))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid field name.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Remove from field map
        // Releases object reference
        fieldDescriptorMap_.erase(fieldName);
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::StringList* FoamX::IApplicationImpl::patchPhysicalTypes()
{
    static const char* functionName =
        "FoamX::IApplicationImpl::patchPhysicalTypes()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate list and return
    return new FoamXServer::StringList
    (
        FoamXWordList
        (
            static_cast<const wordList&>(patchPhysicalTypeDescriptorMap_.toc())
        )
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IApplicationImpl::getPatchPhysicalType
(
    const char* patchPhysicalTypeName,
    IPatchPhysicalTypeDescriptor_out boundaryDescriptor
)
{
    static const char* functionName =
        "FoamX::IApplicationImpl::getPatchPhysicalType"
        "(const char* patchPhysicalTypeName, "
        "IPatchPhysicalTypeDescriptor_out boundaryDescriptor)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // See if we have this object cached
        if (!patchPhysicalTypeDescriptorMap_.found(patchPhysicalTypeName))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid boundary type name "
              + word(patchPhysicalTypeName) + ".",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Return a reference to object
        boundaryDescriptor =
            patchPhysicalTypeDescriptorMap_[patchPhysicalTypeName]->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void  FoamX::IApplicationImpl::findPatchPhysicalType
(
    const char* patchPhysicalTypeName,
    IPatchPhysicalTypeDescriptor_out boundaryDescriptor
)
{
    static const char* functionName =
        "FoamX::IApplicationImpl::findPatchPhysicalType"
        "(const char* patchPhysicalTypeName, "
        "IPatchPhysicalTypeDescriptor_out boundaryDescriptor)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Find the boundary type descriptors with the specified name
        if (patchPhysicalTypeDescriptorMap_.found(patchPhysicalTypeName))
        {
            boundaryDescriptor = 
                patchPhysicalTypeDescriptorMap_.lookup(patchPhysicalTypeName)
                ->_this();
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IApplicationImpl::addPatchPhysicalType
(
    const char* patchPhysicalTypeName,
    IPatchPhysicalTypeDescriptor_out boundaryDescriptor
)
{
    static const char* functionName =
        "FoamX::IApplicationImpl::addPatchPhysicalType"
        "(const char* patchPhysicalTypeName, "
        "IPatchPhysicalTypeDescriptor_out boundaryDescriptor)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {

        // See if we already have a boundary type with this name
        if (patchPhysicalTypeDescriptorMap_.found(patchPhysicalTypeName))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid boundary type name.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Create and initialise a new PatchPhysicalTypeDescriptor object
        IPatchPhysicalTypeDescriptorImpl* pBoundaryDescriptor =
            new IPatchPhysicalTypeDescriptorImpl(patchPhysicalTypeName);

        if (pBoundaryDescriptor == NULL)
        {
            throw FoamXError
            (
                E_FAIL,
                "Couldn't create PatchPhysicalTypeDescriptor object.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Add to map
        patchPhysicalTypeDescriptorMap_.append
        (
            patchPhysicalTypeName,
            pBoundaryDescriptor
        );

        // Return a reference to the PatchPhysicalTypeDescriptor object
        boundaryDescriptor = pBoundaryDescriptor->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IApplicationImpl::deletePatchPhysicalType
(
    const char* patchPhysicalTypeName
)
{
    static const char* functionName =
        "FoamX::IApplicationImpl::deletePatchPhysicalType"
        "(const char* patchPhysicalTypeName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        log << "delete " << patchPhysicalTypeName << endl;

        // See if we have a boundary type with this name
        if (!patchPhysicalTypeDescriptorMap_.found(patchPhysicalTypeName))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid boundary type name.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Remove from map
        // Releases object reference
        patchPhysicalTypeDescriptorMap_.erase(patchPhysicalTypeName);
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::StringList* FoamX::IApplicationImpl::dictionaries()
{
    static const char* functionName =
        "FoamX::IApplicationImpl::dictionaries()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate list and return
    return new FoamXServer::StringList
    (
        FoamXWordList
        (
            static_cast<const wordList&>(dictTypeDescriptorMap_.toc())
        )
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IApplicationImpl::getDictionary
(
    const char* dictName,
    ITypeDescriptor_out dictTypeDescriptor
)
{
    static const char* functionName =
        "FoamX::IApplicationImpl::getDictionary"
        "(const char* dictName, ITypeDescriptor_out dictTypeDescriptor)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {

        // See if we have this dictionary's type descriptor cached
        if (!dictTypeDescriptorMap_.found(dictName))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid dictionary name '" + word(dictName) + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Return a reference to the TypeDescriptor object
        dictTypeDescriptor = dictTypeDescriptorMap_[dictName]->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IApplicationImpl::addDictionary
(
    const char* dictName,
    ITypeDescriptor_out dictTypeDescriptor
)
{
    static const char* functionName =
        "FoamX::IApplicationImpl::addDictionary"
        "(const char* dictName, ITypeDescriptor_out dictTypeDescriptor)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {

        // See if we already have a dictionary with this name
        if (dictTypeDescriptorMap_.found(dictName))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid dictionary name.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Create and initialise a new TypeDescriptor object for a new
        // dictionary
        ITypeDescriptorImpl* pTypeDescriptor = new ITypeDescriptorImpl
        (
            dictName,
            Type_Dictionary,
            dictName
        );

        if (pTypeDescriptor == NULL)
        {
            throw FoamXError
            (
                E_FAIL,
                "Couldn't create Dictionary TypeDescriptor object.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Add to dictionary map
        dictTypeDescriptorMap_.append(dictName, pTypeDescriptor);

        // Return a reference to the TypeDescriptor object
        dictTypeDescriptor = pTypeDescriptor->_this();
    }
    CATCH_ALL(functionName);

    return;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IApplicationImpl::deleteDictionary(const char* dictName)
{
    static const char* functionName =
        "FoamX::IApplicationImpl::deleteDictionary(const char* dictName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {

        // See if we have a dictionary with this name
        if (!dictTypeDescriptorMap_.found(dictName))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid dictionary name.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Remove from dictionary map
        // Releases object reference
        dictTypeDescriptorMap_.erase(dictName);
    }
    CATCH_ALL(functionName);
}


// ************************************************************************* //
