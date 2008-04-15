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
#include "IGeometricFieldDescriptorImpl.H"
#include "ITypeDescriptorImpl.H"
#include "IGeometryDescriptorImpl.H"
#include "DimensionSet.H"
#include "LogEntry.H"
#include "Paths.H"

// Namespaces
#include "FoamXNameSpaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IGeometricFieldDescriptorImpl::IGeometricFieldDescriptorImpl
(
    const char* fieldName
)
:
    fieldName_(fieldName),
    fieldDescription_(fieldName),
    typeDescriptor_(NULL),
    geometryDescriptor_(NULL)
{
    static const char* functionName =
        "FoamX::IGeometricFieldDescriptorImpl::IGeometricFieldDescriptorImpl"
        "(const char* fieldName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Initialise dimension set.
    dimensionSet_ == dimless;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IGeometricFieldDescriptorImpl::IGeometricFieldDescriptorImpl
(
    const word& fieldName,
    const dictionary& fieldDict,
    const ObjRefHashTable<ITypeDescriptorImpl*>& foamTypes,
    const ObjRefHashTable<IGeometryDescriptorImpl*>& geometryTypes
)
:
    fieldName_(fieldName),
    typeDescriptor_(NULL),
    geometryDescriptor_(NULL)
{
    static const char* functionName =
        "FoamX::IGeometricFieldDescriptorImpl::IGeometricFieldDescriptorImpl"
        "(const word& fieldName, const dictionary& fieldDict)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Extract the field type name and description.
        fieldDict.lookup("description") >> fieldDescription_;
        fieldDict.lookup("fieldType") >> fieldTypeName_;
        fieldDict.lookup("geometryType") >> geometryTypeName_;

        // Extract the dimension set.
        dimensionSet_ == dimensionSet(fieldDict.lookup("dimensions"));

        // Get the field descriptor for the specified field type.
        if (!foamTypes.found(fieldTypeName_))
        {
            throw FoamXError
            (
                E_FAIL,
                "Type descriptor not found for "
              + fieldTypeName_,
                functionName,
                __FILE__, __LINE__
            );
        }
        typeDescriptor_ = foamTypes.find(fieldTypeName_)();


        // Get the field descriptor for the specified field type.
        word fieldName = fieldTypeName_ + "Field";
        if (!foamTypes.found(fieldName))
        {
            throw FoamXError
            (
                E_FAIL,
                "Type descriptor not found for field type "
              + fieldName,
                functionName,
                __FILE__, __LINE__
            );
        }
        fieldTypeDescriptor_ = foamTypes.find(fieldName)();


        // Get the geometry descriptor for the specified geometry type.
        if (!geometryTypes.found(geometryTypeName_))
        {
            throw FoamXError
            (
                E_FAIL,
                "Geometry descriptor not found for specified geometry type "
              + geometryTypeName_,
                functionName,
                __FILE__, __LINE__
            );
        }
        geometryDescriptor_ = geometryTypes.find(geometryTypeName_)();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IGeometricFieldDescriptorImpl::~IGeometricFieldDescriptorImpl()
{
    static const char* functionName =
        "FoamX::IGeometricFieldDescriptorImpl::"
        "~IGeometricFieldDescriptorImpl()";

    LogEntry log(functionName, __FILE__, __LINE__);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::ITypeDescriptor_ptr 
FoamX::IGeometricFieldDescriptorImpl::typeDescriptor()
{
    static const char* functionName =
        "FoamX::IGeometricFieldDescriptorImpl::typeDescriptor()";

    LogEntry log(functionName, __FILE__, __LINE__);

    ITypeDescriptor_ptr pTypeDesc = ITypeDescriptor::_nil();

    if (typeDescriptor_ != NULL)
    {
        // Return a duplictate reference to the TypeDescriptor object.
        pTypeDesc = typeDescriptor_->_this();
    }

    return pTypeDesc;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::ITypeDescriptor_ptr 
FoamX::IGeometricFieldDescriptorImpl::fieldTypeDescriptor()
{
    static const char* functionName =
        "FoamX::IGeometricFieldDescriptorImpl::fieldTypeDescriptor()";

    LogEntry log(functionName, __FILE__, __LINE__);

    ITypeDescriptor_ptr pTypeDesc = ITypeDescriptor::_nil();

    if (typeDescriptor_ != NULL)
    {
        // Return a duplictate reference to the TypeDescriptor object.
        pTypeDesc = fieldTypeDescriptor_->_this();
    }

    return pTypeDesc;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IGeometryDescriptor_ptr
FoamX::IGeometricFieldDescriptorImpl::geometryDescriptor()
{
    static const char* functionName =
        "FoamX::IGeometricFieldDescriptorImpl::geometryDescriptor()";

    LogEntry log(functionName, __FILE__, __LINE__);

    IGeometryDescriptor_ptr pGeomDesc = IGeometryDescriptor::_nil();

    if (geometryDescriptor_ != NULL)
    {
        // Return a duplictate reference to the GeometryDescriptor object.
        pGeomDesc = geometryDescriptor_->_this();
    }

    return pGeomDesc;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::IGeometricFieldDescriptorImpl::name()
{
    static const char* functionName =
        "FoamX::IGeometricFieldDescriptorImpl::name()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate string and return.
    return CORBA::string_dup(fieldName_.c_str());
}

void FoamX::IGeometricFieldDescriptorImpl::name(const char* newName)
{
    static const char* functionName =
        "FoamX::IGeometricFieldDescriptorImpl::name(const char* newName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    fieldName_ = newName;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::IGeometricFieldDescriptorImpl::description()
{
    static const char* functionName =
        "FoamX::IGeometricFieldDescriptorImpl::description()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate string and return.
    return CORBA::string_dup(fieldDescription_.c_str());
}

void FoamX::IGeometricFieldDescriptorImpl::description
(
    const char* newDescription
)
{
    static const char* functionName =
        "FoamX::IGeometricFieldDescriptorImpl::"
        "description(const char* newDescription)";

    LogEntry log(functionName, __FILE__, __LINE__);

    fieldDescription_ = newDescription;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::IGeometricFieldDescriptorImpl::fieldTypeName()
{
    static const char* functionName =
        "FoamX::IGeometricFieldDescriptorImpl::fieldTypeName()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate string and return.
    return CORBA::string_dup(fieldTypeName_.c_str());
}

void FoamX::IGeometricFieldDescriptorImpl::fieldTypeName(const char* newName)
{
    static const char* functionName =
        "FoamX::IGeometricFieldDescriptorImpl::"
        "fieldTypeName(const char* newName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    fieldTypeName_ = newName;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::IGeometricFieldDescriptorImpl::geometryTypeName()
{
    static const char* functionName =
        "FoamX::IGeometricFieldDescriptorImpl::geometryTypeName()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate string and return.
    return CORBA::string_dup(geometryTypeName_.c_str());
}

void FoamX::IGeometricFieldDescriptorImpl::geometryTypeName(const char* newName)
{
    static const char* functionName =
        "FoamX::IGeometricFieldDescriptorImpl::"
        "geometryTypeName(const char* newName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    geometryTypeName_ = newName;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::DimensionSet FoamX::IGeometricFieldDescriptorImpl::dimensions()
{
    static const char* functionName =
        "FoamX::IGeometricFieldDescriptorImpl::dimensions()";

    LogEntry log(functionName, __FILE__, __LINE__);

    return dimensionSet_;
}

void FoamX::IGeometricFieldDescriptorImpl::dimensions
(
    const DimensionSet& newDimSet
)
{
    static const char* functionName =
        "FoamX::IGeometricFieldDescriptorImpl::dimensions"
        "(const DimensionSet& newDimSet)";

    LogEntry log(functionName, __FILE__, __LINE__);

    dimensionSet_ = newDimSet;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IGeometricFieldDescriptorImpl::save(DictionaryWriter& dict)
{
    static const char* functionName =
        "FoamX::IGeometricFieldDescriptorImpl::save(DictionaryWriter& dict)";

    LogEntry log(functionName, __FILE__, __LINE__);


    // Save this field information to the dictionary.
    dict.startSubDict(fieldName_);

    dict.writeEntry("description", fieldDescription_);
    dict.writeEntry("fieldType", fieldTypeName_);
    dict.writeEntry("geometryType", geometryTypeName_);
    dict.writeEntry("dimensions", dimensionSet_);

    dict.endSubDict();
}


// ************************************************************************* //
