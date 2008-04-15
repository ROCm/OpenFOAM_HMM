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
#include "OSspecific.H"

// FoamX header files.
#include "FoamX.H"
#include "FoamXErrors.H"
#include "IGeometricFieldImpl.H"
#include "IDictionaryEntryImpl.H"
#include "LogEntry.H"
#include "Paths.H"

// Namespaces
#include "FoamXNameSpaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IGeometricFieldImpl::IGeometricFieldImpl
(
    IGeometricFieldDescriptor_ptr fieldDescriptor,
    IFoamProperties_ptr foamProperties
)
:
    internalFieldValue_(NULL)
{
    static const char* functionName =
        "FoamX::IGeometricFieldImpl::IGeometricFieldImpl"
        "(IGeometricFieldDescriptor_ptr fieldDescriptor, "
        "IFoamProperties_ptr foamProperties)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Take references to the FieldDescriptor and FoamProperties objects.
        if (CORBA::is_nil(foamProperties))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid FoamProperties reference.",
                functionName,
                __FILE__, __LINE__
            );
        }
        if (CORBA::is_nil(fieldDescriptor))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid FieldDescriptor reference.",
                functionName,
                __FILE__, __LINE__
            );
        }
        foamProperties_ = IFoamProperties::_duplicate(foamProperties);
        fieldDescriptor_ = IGeometricFieldDescriptor::_duplicate
        (
            fieldDescriptor
        );

        // Get the field name.
        fieldName_ = fieldDescriptor_->name();

        // Create the internal field value.
        ITypeDescriptor_var fieldType = fieldDescriptor_->fieldTypeDescriptor();
        internalFieldValue_ = new IDictionaryEntryImpl(fieldType);
        if (internalFieldValue_ == NULL)
        {
            throw FoamXError
            (
                E_FAIL,
                "Failed to create internal field value "
                "dictionary entry object.",
                functionName,
                __FILE__, __LINE__
            );
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IGeometricFieldImpl::~IGeometricFieldImpl()
{
    static const char* functionName =
        "FoamX::IGeometricFieldImpl::~IGeometricFieldImpl()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Release internal value and reference level value dictionary entry objects
    if (internalFieldValue_  != NULL)
    {
        internalFieldValue_->_remove_ref();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::IGeometricFieldImpl::name()
{
    static const char* functionName =
        "FoamX::IGeometricFieldImpl::name()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Get name from fieldDescriptor object.
    if (!CORBA::is_nil(fieldDescriptor_))
    {
        return fieldDescriptor_->name();
    }
    else
    {
        return CORBA::string_dup("");
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IGeometricFieldImpl::getInternalFieldValue
(
    IDictionaryEntry_out internalFieldValue
)
{
    static const char* functionName =
        "FoamX::IGeometricFieldImpl::getInternalFieldValue"
        "(IDictionaryEntry_out internalFieldValue)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        if (internalFieldValue_ == NULL)
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Invalid internal field value object.",
                functionName,
                __FILE__, __LINE__
            );
        }

        internalFieldValue = internalFieldValue_->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IGeometricFieldImpl::getPatchFieldParameters
(
    const char* patchName,
    IDictionaryEntry_out patchFields
)
{
    static const char* functionName =
        "FoamX::IGeometricFieldImpl::getPatchFieldParameters"
        "(const char* patchName, IDictionaryEntry_out patchFields)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // See if this patch requires extra parameters.
        if
        (
            patchFields_.found(patchName)
         && patchFields_[patchName]->subElementPtrs().size() > 1
        )
        {
            patchFields = patchFields_[patchName]->_this();
        }
        else
        {
            patchFields = IDictionaryEntry::_nil();
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Set the patch field type for the specified patch. Note that the patch
// field type parameter is the true Foam name and not the patch
// field descriptor key.

void FoamX::IGeometricFieldImpl::addPatch
(
    const char* patchName,
    const char* patchFieldType
)
{
    static const char* functionName =
        "FoamX::IGeometricFieldImpl::addPatch"
        "(const char* patchName, const char* patchFieldType)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        log << "Adding patch " << patchName
            << " to field " << fieldName_ << endl;

        // Check patch name.
        if (patchFieldNames_.found(patchName))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid patch name.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Lookup this type and see if additional parameters are required.
        ITypeDescriptor_var pfDesc;
        foamProperties_->findPatchFieldType(patchFieldType, pfDesc.out());

        // Make sure that the patch field type name is valid.
        if (CORBA::is_nil(pfDesc))
        {
            throw FoamXError
            (
                E_FAIL,
                "Invalid patch field type name '" + word(patchFieldType) + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Set patch field type.
        patchFieldNames_.insert(patchName, patchFieldType);

        // Get the type descriptor for the extra parameters. May be nil.
        if (!CORBA::is_nil(pfDesc))
        {
            // Create a default DictionaryEntry object.
            IDictionaryEntryImpl* patchDictEntryPtr = new IDictionaryEntryImpl
            (
                pfDesc
            );

            if (patchDictEntryPtr == NULL)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Failed to create IDictionaryEntryImpl object.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Bind the field type descriptor to the parameter
            // dictionary entry object.
            patchDictEntryPtr->bindFieldType
            (
                fieldDescriptor_->fieldTypeDescriptor()
            );

            // Add to map.
            patchFields_.insert(patchName, patchDictEntryPtr);
        }

        log << "Added patch " << patchName
            << " to field " << fieldName_ << endl;
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IGeometricFieldImpl::deletePatch(const char* patchName)
{
    static const char* functionName =
        "FoamX::IGeometricFieldImpl::deletePatch(const char* patchName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Check patch name.
        if (!patchFieldNames_.found(patchName))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid patch name.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Remove from patch field type map.
        patchFieldNames_.erase(patchFieldNames_.find(patchName));

        // Remove from patch field parameter map if necessary.
        if (patchFields_.found(patchName))
        {
            HashTable<IDictionaryEntryImpl*>::iterator iter =
                patchFields_.find(patchName);
            patchFields_.erase(iter);     // Releases reference.
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IGeometricFieldImpl::setPatchFieldType
(
    const char* patchName,
    const char* patchFieldType
)
{
    static const char* functionName =
        "FoamX::IGeometricFieldImpl::setPatchFieldType"
        "(const char* patchName, const char* patchFieldType)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Set patch field type if it has changed.
        if
        (
            !patchFieldNames_.found(patchName)
          || patchFieldType != patchFieldNames_.find(patchName)()
        )
        {
            deletePatch(patchName);
            addPatch(patchName, patchFieldType);
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

CORBA::Boolean FoamX::IGeometricFieldImpl::modified()
{
    static const char* functionName =
        "FoamX::ICaseServerImpl::modified()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        if (internalFieldValue_->modified()) return true;

        for
        (
            ObjRefHashTable<IDictionaryEntryImpl*>::iterator iter =
                patchFields_.begin();
            iter != patchFields_.end();
            ++iter
        )
        {
            if (iter()->modified()) return true;
        }
    }
    CATCH_ALL(functionName);

    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IGeometricFieldImpl::load(const dictionary& fieldDict)
{
    static const char* functionName =
        "FoamX::IGeometricFieldImpl::load(const dictionary& fieldDict)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        log << "loading field " << fieldName_ << endl;

        log << "Reading " << fieldDict.name()
            << " start line " << fieldDict.startLineNumber()
            << " end line " << fieldDict.endLineNumber() << endl;

        // Make sure the field dictionary contains the required entries.
        if
        (
            !fieldDict.found("internalField")
         || !fieldDict.found("boundaryField")
        )
        {
            throw FoamXError
            (
                E_FAIL,
                "Malformed field dictionary.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Read the internal field value.
        internalFieldValue_->load(fieldDict.lookup("internalField"));

        // Open boundary field dictionary.
        const dictionary& boundaryFieldDict(fieldDict.subDict("boundaryField"));

        // Loop over all defined patch names and read the patch field values.
        wordList patchNames = patchFieldNames_.toc();
        forAll(patchNames, i)
        {
            word patchName(patchNames[i]);

            if (!boundaryFieldDict.found(patchName))
            {
                // Patch dictionary not found. Bugger.
                log << "Warning : Patch dictionary '" << patchName
                    << "' not found in field dictionary '" << fieldName_
                    << "'." << endl;

                WarningIn(functionName)
                    << "Patch dictionary '" << patchName
                    << "' not found in field dictionary '" << fieldName_
                    << "'." << endl;
                continue;
            }

            // Open patch dictionary.
            const dictionary& patchDict(boundaryFieldDict.subDict(patchName));

            // Make sure the patch dictionary contains the required entries.
            if (!patchDict.found("type"))
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Patch field type not specified for patch '" + patchName
                  + "' in field dictionary '" + fieldName_ + "'.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Read the patch field type. This is the actual type not the key.
            word patchFieldType(patchDict.lookup("type"));

            // Get the PatchFieldDescriptor object. Search by type.
            ITypeDescriptor_var pfDesc;
            foamProperties_->findPatchFieldType
            (
                patchFieldType.c_str(),
                pfDesc.out()
            );

            // Make sure that the patch field type name is valid.
            if (CORBA::is_nil(pfDesc))
            {
                log << "Warning : Invalid patch field type '" << patchFieldType
                    << "' for patch '" << patchName
                    << "' in field dictionary '" << fieldName_ << "'."
                    << endl;

                WarningIn(functionName)
                    << "Invalid patch field type '" << patchFieldType
                    << "' for patch '" << patchName
                    << "' in field dictionary '" << fieldName_ << "'."
                    << endl;

                continue;
            }

            // Check that the patch field type is valid for the boundary type
            // specified for this patch.
            word boundaryPatchFieldType = patchFieldNames_[patchName];
            if (patchFieldType != boundaryPatchFieldType)
            {
                log << "Warning : Incorrect patch field type '"
                    << patchFieldType << "' for patch '" << patchName << "'."
                    << endl << "          Boundary condition specifies '"
                    << boundaryPatchFieldType
                    << "' for field '" << fieldName_ << "'."
                    << endl;

                WarningIn(functionName)
                    << "Incorrect patch field type '"
                    << patchFieldType << "' for patch '" << patchName << "'."
                    << endl << "          Boundary condition specifies '"
                    << boundaryPatchFieldType
                    << "' for field '" << fieldName_ << "'."
                    << endl;
            }

            // Set the patch field type. Poss over-ride default?
            patchFieldNames_[patchName] = patchFieldType;

            // Release any previous patch field parameter object.
            if (patchFields_.found(patchName))
            {
                HashTable<IDictionaryEntryImpl*>::iterator iter =
                    patchFields_.find(patchName);
                patchFields_.erase(iter);
            }

            // Get the type descriptor for the extra parameters. May be nil.
            if (!CORBA::is_nil(pfDesc))
            {
                // Create and initialise a new IDictionaryEntryImpl
                // object for the parameters.
                IDictionaryEntryImpl* patchDictEntryPtr =
                new IDictionaryEntryImpl
                (
                    pfDesc
                );
                if (patchDictEntryPtr == NULL)
                {
                    throw FoamXError
                    (
                        E_FAIL,
                        "Failed to create patch field parameter "
                        "DictionaryEntry Object.",
                        functionName,
                        __FILE__, __LINE__
                    );
                }

                // Bind the field type descriptor to the parameter
                // dictionary entry object.
                patchDictEntryPtr->bindFieldType
                (
                    fieldDescriptor_->fieldTypeDescriptor()
                );

                // Load values from dictionary.
                // The extra parameters are defined as dictionary
                // entries. Pass the entire patch dictionary through.
                // Tolerate non-existent dictionary entries.
                // Any missing parameters will get a default value.

                patchDictEntryPtr->load
                (
                    boundaryFieldDict.subDict(patchName),
                    true
                );

                // Add to map.
                patchFields_.insert(patchName, patchDictEntryPtr);
            }
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IGeometricFieldImpl::save
(
    DictionaryWriter& dictWriter,
    const wordList& patchNames
)
{
    static const char* functionName =
        "FoamX::IGeometricFieldImpl::save"
        "(DictionaryWriter& dictWriter, const wordList& patchNames)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Validate the dictionary entry items.
        internalFieldValue_->validate();

        // Construct class name (eg, "volScalarField").
        word className
        (
            word(fieldDescriptor_->geometryDescriptor()->name())
          & word(fieldDescriptor_->fieldTypeDescriptor()->name())
        );

        dictWriter.writeHeader
        (
            "Field Dictionary",
            className
        );

        dictWriter.writeEndl();

        // Get field dimensions from field descriptor.
        DimensionSet fieldDimension = fieldDescriptor_->dimensions();
        dictWriter.writeEntry("dimensions", fieldDimension);
        dictWriter.writeEndl();

        // Write internal field and reference values.
        dictWriter.writeKeyword("internalField");
        internalFieldValue_->save(dictWriter, false);
        dictWriter.endEntry();
        dictWriter.writeEndl();

        // Open boundary field sub-dictionary.
        dictWriter.startSubDict("boundaryField");

        // Loop over all patchs.
        forAll(patchNames, i)
        {
            dictWriter.writeKeyword(patchNames[i]);

            // Write value only.
            patchFields_[patchNames[i]]->save(dictWriter, false);
            dictWriter.writeEndl();

            if (i < patchNames.size()-1)
            {
                dictWriter.writeEndl();
            }
        }

        dictWriter.endSubDict();
        dictWriter.writeEndl();
        dictWriter.writeEndBar();
    }
    catch (FoamXError& ex)
    {
        // Save failed. Remove the duff file.
        if (exists(dictWriter.pathName()))
        {
            rm(dictWriter.pathName());
        }

        // Bounce exception up to client.
        throw ex;
    }
    catch (...)
    {
        // Save failed. Remove the duff file.
        if (exists(dictWriter.pathName()))
        {
            rm(dictWriter.pathName());
        }

        throw FoamXError
        (
            E_UNEXPECTED,
            "Unexpected error.",
            functionName,
            __FILE__, __LINE__
        );
    }
}


// ************************************************************************* //
