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

// Foam header files.
#include "SLList.H"
#include "OSspecific.H"
#include "IFstream.H"

// Project header files.
#include "ITypeDescriptorImpl.H"
#include "FoamX.H"
#include "FoamXErrors.H"
#include "FoamXAnyList.H"
#include "LogEntry.H"

// Namespaces
#include "FoamXNameSpaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::ITypeDescriptorImpl::ITypeDescriptorImpl
(
    const word& name,
    const FoamXType& type,
    const string& parentPath
)
:
    type_(type),
    name_(name),
    path_(parentPath + ':' + name),
    optional_(false),
    visible_(true),
    editable_(true),
    numElements_(0),
    default_(NULL)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::ITypeDescriptorImpl"
        "(const word& name, const FoamXType& type, "
        "const string& parentPath)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        log << "Constructing TypeDescriptor for " << path_ << "." << endl;

        if (FoamXTypes::isNumber(type_))
        {
            minValue_.setType(type_);
            minValue_.setMin();

            maxValue_.setType(type_);
            maxValue_.setMax();
        }
    }
    CATCH_ALL(functionName);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::ITypeDescriptorImpl::ITypeDescriptorImpl
(
    const word& name,
    const string& parentPath,
    const dictionary& typeDict,
    const dictionary& foamTypesDict
)
:
    type_(Type_Undefined),
    name_(name),
    path_(parentPath + ':' + name),
    optional_(false),
    visible_(true),
    editable_(true),
    numElements_(0),
    default_(NULL)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::ITypeDescriptorImpl"
        "(const word& name, const string& parentPath, "
        "const dictionary& typeDict, "
        "const dictionary& foamTypesDict)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        log << "Constructing TypeDescriptor for " << path_ << "." << endl;

        log << "Reading " << typeDict.name()
            << " start line " << typeDict.startLineNumber()
            << " end line " << typeDict.endLineNumber() << endl;

        // Load the type descriptor information from the dictionary.
        load(typeDict, foamTypesDict);
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::ITypeDescriptorImpl::ITypeDescriptorImpl
(
    const word& name,
    const string& parentPath,
    const entry& typeEntry,
    const dictionary& foamTypesDict
)
:
    type_(Type_Undefined),
    name_(name),
    path_(parentPath + ':' + name),
    optional_(false),
    visible_(true),
    editable_(true),
    numElements_(0),
    default_(NULL)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::ITypeDescriptorImpl"
        "(const word& name, const string& parentPath, const entry& typeEntry, "
        "const dictionary& foamTypesDict)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        log << "Constructing TypeDescriptor for " << path_ << "." << endl;

        if (typeEntry.isDict())
        {
            // Load the type descriptor information from the dictionary.
            load(typeEntry.dict(), foamTypesDict);
        }
        else if (exists(parentPath/name_ + ".cfg"))
        {
            fileName dictPathName = parentPath/name_ + ".cfg";
        
            dictionary typeDict((IFstream(dictPathName)()));

            log << "Reading dictionary " << typeDict.name()
                << " start line " << typeDict.startLineNumber()
                << " end line " << typeDict.endLineNumber() << endl;

            // Load the type descriptor information from the dictionary.
            load(typeDict, foamTypesDict);
        }
        else if (foamTypesDict.found(name_))
        {
            // Load the type descriptor information from the dictionary.
            load(foamTypesDict.subDict(name_), foamTypesDict);
        }
        else
        {
            throw FoamXError
            (
                E_FAIL,
                "Definition for type '" + name_ + "', " + path_ + " could not be found in "
              + typeEntry.name() + " or in foamTypesDict",
                functionName,
                __FILE__, __LINE__
            );
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::ITypeDescriptorImpl::ITypeDescriptorImpl
(
    ITypeDescriptor_ptr typeDesc
)
:
    type_(Type_Undefined),
    optional_(false),
    visible_(true),
    editable_(true),
    numElements_(0),
    default_(NULL)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::ITypeDescriptorImpl"
        "(ITypeDescriptor_ptr typeDesc)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        log << "Constructing TypeDescriptor." << endl;

        // Set basic type.
        type_ = typeDesc->type();

        // Copy common data
        name_ = typeDesc->name();
        path_ = typeDesc->path();
        displayName_ = typeDesc->displayName();
        description_ = typeDesc->description();
        comment_ = typeDesc->comment();
        category_ = typeDesc->category();
        helpURL_ = typeDesc->helpURL();
        iconURL_ = typeDesc->iconURL();
        optional_ = typeDesc->optional();
        visible_  = typeDesc->visible();
        editable_ = typeDesc->editable();

        if (!isCompoundType())
        {
            if (FoamXTypes::isNumber(type_))
            {
                // Copy non-compound specific information.
                FoamXServer::FoamXAny_var min = typeDesc->minValue();
                minValue_.setType(type_);
                minValue_.setValue(min);

                FoamXServer::FoamXAny_var max = typeDesc->maxValue();
                maxValue_.setType(type_);
                maxValue_.setValue(max);
            }

            lookupDict_ = typeDesc->lookupDict();

            FoamXServer::FoamXAnyList_var anyList = typeDesc->valueList();

            valueList_.setSize(anyList->length());

            for (label i=0; i<valueList_.size(); i++)
            {
                valueList_[i].setType(type_);
                valueList_[i].setValue(anyList[ulong(i)]);
            }
        }
        else
        {
            // Copy compound specific information.
            dictionaryPath_ = typeDesc->dictionaryPath();

            numElements_ = typeDesc->numElements();
            elementLabels_ = StringList_var(typeDesc->elementLabels());

            // Get list of the sub types from the type descriptor.
            // Auto release.
            TypeDescriptorList_var pSubTypes = typeDesc->subTypes();

            // Loop over all sub types.
            for (unsigned int i = 0; i <pSubTypes->length(); i++)
            {
                // Get the TypeDescriptor for this entry.
                ITypeDescriptor_var pSubTypeDesc =
                    FoamXServer::ITypeDescriptor::_duplicate(pSubTypes[i]);

                // Construct the type descriptor.
                ITypeDescriptorImpl* pTypeDescriptor =
                    new ITypeDescriptorImpl(pSubTypeDesc);

                if (pTypeDescriptor == NULL)
                {
                    throw FoamXError
                    (
                        E_FAIL,
                        "Couldn't create duplicate type descriptor for "
                      + word(pSubTypeDesc->name()),
                        functionName,
                        __FILE__, __LINE__
                    );
                }

                // Add to map and list.
                word name(pSubTypeDesc->name());
                subTypes_.append(pTypeDescriptor);
            }
        }

        if (typeDesc->hasDefaultValue())
        {
            // Construct default_ from the given type descriptor
            // which automatically sets the values to those contained in the
            // default_ entry of the type descriptor
            default_ = new IDictionaryEntryImpl(typeDesc);
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::ITypeDescriptorImpl::~ITypeDescriptorImpl()
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::~ITypeDescriptorImpl()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Release default dictionary entry object
    if (default_ != NULL)
    {
        default_->_remove_ref();
    }

    // Note : Sub-type objects are released in ObjRefHashTable destructor.
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

CORBA::Boolean FoamX::ITypeDescriptorImpl::isPrimitiveType()
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::isPrimitiveType()";

    LogEntry log(functionName, __FILE__, __LINE__);

    return FoamXTypes::isPrimitive(type_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

CORBA::Boolean FoamX::ITypeDescriptorImpl::isCompoundType()
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::isCompoundType()";

    LogEntry log(functionName, __FILE__, __LINE__);

    return FoamXTypes::isCompound(type_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::FoamXType FoamX::ITypeDescriptorImpl::type()
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::type()";

    LogEntry log(functionName, __FILE__, __LINE__);

    return type_;
}

void FoamX::ITypeDescriptorImpl::type(FoamXType type)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::type(FoamXType type)";

    LogEntry log(functionName, __FILE__, __LINE__);

    type_ = type;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::ITypeDescriptorImpl::path()
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::path()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate and return.
    return CORBA::string_dup(path_.c_str());
}

void FoamX::ITypeDescriptorImpl::path(const char* path)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::path(const char* path)";

    LogEntry log(functionName, __FILE__, __LINE__);

    path_ = path;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::ITypeDescriptorImpl::name()
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::name()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate and return.
    return CORBA::string_dup(name_.c_str());
}

void FoamX::ITypeDescriptorImpl::name(const char* name)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::name(const char* name)";

    LogEntry log(functionName, __FILE__, __LINE__);

    name_ = name;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::ITypeDescriptorImpl::displayName()
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::displayName()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate and return.
    return CORBA::string_dup(displayName_.c_str());
}

void FoamX::ITypeDescriptorImpl::displayName(const char* name)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::displayName(const char* name)";

    LogEntry log(functionName, __FILE__, __LINE__);

    displayName_ = name;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::ITypeDescriptorImpl::description()
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::description()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate and return.
    return CORBA::string_dup(description_.c_str());
}

void FoamX::ITypeDescriptorImpl::description(const char* desc)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::description(const char* desc)";

    LogEntry log(functionName, __FILE__, __LINE__);

    description_ = desc;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::ITypeDescriptorImpl::category()
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::category()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate and return.
    return CORBA::string_dup(category_.c_str());
}

void FoamX::ITypeDescriptorImpl::category(const char* cat)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::category(const char* cat)";

    LogEntry log(functionName, __FILE__, __LINE__);

    category_ = cat;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::ITypeDescriptorImpl::comment()
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::comment()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate and return.
    return CORBA::string_dup(comment_.c_str());
}

void FoamX::ITypeDescriptorImpl::comment(const char* com)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::comment(const char* com)";

    LogEntry log(functionName, __FILE__, __LINE__);

    comment_ = com;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::ITypeDescriptorImpl::helpURL()
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::helpURL()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate and return.
    return CORBA::string_dup(helpURL_.c_str());
}

void FoamX::ITypeDescriptorImpl::helpURL(const char* url)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::helpURL(const char* url)";

    LogEntry log(functionName, __FILE__, __LINE__);

    helpURL_ = url;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::ITypeDescriptorImpl::iconURL()
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::iconURL()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate and return.
    return CORBA::string_dup(iconURL_.c_str());
}

void FoamX::ITypeDescriptorImpl::iconURL(const char* url)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::iconURL(const char* url)";

    LogEntry log(functionName, __FILE__, __LINE__);

    iconURL_ = url;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

CORBA::Boolean FoamX::ITypeDescriptorImpl::optional()
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::optional()";

    LogEntry log(functionName, __FILE__, __LINE__);

    return optional_;
}

void FoamX::ITypeDescriptorImpl::optional(CORBA::Boolean opt)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::optional(CORBA::Boolean opt)";

    LogEntry log(functionName, __FILE__, __LINE__);

    optional_ = opt;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

CORBA::Boolean FoamX::ITypeDescriptorImpl::visible()
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::visible()";

    LogEntry log(functionName, __FILE__, __LINE__);

    return visible_;
}

void FoamX::ITypeDescriptorImpl::visible(CORBA::Boolean vis)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::visible(CORBA::Boolean vis)";

    LogEntry log(functionName, __FILE__, __LINE__);

    visible_ = vis;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

CORBA::Boolean FoamX::ITypeDescriptorImpl::editable()
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::editable()";

    LogEntry log(functionName, __FILE__, __LINE__);

    return editable_;
}

void FoamX::ITypeDescriptorImpl::editable(CORBA::Boolean edit)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::editable(CORBA::Boolean edit)";

    LogEntry log(functionName, __FILE__, __LINE__);

    editable_ = edit;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::FoamXAny* FoamX::ITypeDescriptorImpl::minValue()
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::minValue()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate and return.
    return new FoamXServer::FoamXAny(minValue_);
}

void FoamX::ITypeDescriptorImpl::minValue(const FoamXServer::FoamXAny& newValue)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::minValue(const FoamXServer::Any& newValue)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Set the value of the any object.
        // Will throw an exception if the incoming any is not compatible.
        minValue_.setValue(newValue);
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::FoamXAny* FoamX::ITypeDescriptorImpl::maxValue()
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::maxValue()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate and return.
    return new FoamXServer::FoamXAny(maxValue_);
}

void FoamX::ITypeDescriptorImpl::maxValue(const FoamXServer::FoamXAny& newValue)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::maxValue(const FoamXServer::Any& newValue)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Set the value of the any object.
        // Will throw an exception if the incoming any is not compatible.
        maxValue_.setValue(newValue);
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::FoamXAnyList* FoamX::ITypeDescriptorImpl::valueList()
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::valueList()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Allocate an AnyList and return.
    FoamXServer::FoamXAnyList* pList = new FoamXServer::FoamXAnyList();

    // Set list length.
    pList->length(valueList_.size());

    forAll(valueList_, i)
    {
        (*pList)[i] = valueList_[i];
    }

    return pList;
}

void FoamX::ITypeDescriptorImpl::valueList
(
    const FoamXServer::FoamXAnyList& anyList
)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::valueList"
        "(const FoamXServer::FoamXAnyList& anyList)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Set the value of the any objects.
        // Will throw an exception if the incoming any is not compatible.
        valueList_.setSize(anyList.length());
        forAll(valueList_, i)
        {
            valueList_[i].setType(type_);
            valueList_[i].setValue(anyList[i]);
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::ITypeDescriptorImpl::lookupDict()
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::lookupDict()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate and return.
    return CORBA::string_dup(lookupDict_.c_str());
}

void FoamX::ITypeDescriptorImpl::lookupDict(const char* lookup)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::lookupDict(const char* lookup)";

    LogEntry log(functionName, __FILE__, __LINE__);

    lookupDict_ = lookup;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

char* FoamX::ITypeDescriptorImpl::dictionaryPath()
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::dictionaryPath()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate and return.
    return CORBA::string_dup(dictionaryPath_.c_str());
}

void FoamX::ITypeDescriptorImpl::dictionaryPath(const char* path)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::dictionaryPath(const char* path)";

    LogEntry log(functionName, __FILE__, __LINE__);

    dictionaryPath_ = path;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

CORBA::Long FoamX::ITypeDescriptorImpl::numElements()
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::numElements()";

    LogEntry log(functionName, __FILE__, __LINE__);

    return numElements_;
}

void FoamX::ITypeDescriptorImpl::numElements(CORBA::Long numElements)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::numElements(CORBA::Long numElements)";

    LogEntry log(functionName, __FILE__, __LINE__);

    numElements_ = numElements;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::StringList* FoamX::ITypeDescriptorImpl::elementLabels()
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::elementLabels()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate list and return.
    return new FoamXServer::StringList(elementLabels_);
}

void FoamX::ITypeDescriptorImpl::elementLabels
(
    const FoamXServer::StringList& newList
)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::elementLabels"
        "(const FoamXServer::StringList& newList)";

    LogEntry log(functionName, __FILE__, __LINE__);

    elementLabels_ = newList;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::TypeDescriptorList* FoamX::ITypeDescriptorImpl::subTypes()
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::subTypes()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Construct a duplicate TypeDescriptorList.
    TypeDescriptorList* typeList = new TypeDescriptorList();
    typeList->length(subTypes_.size());
    int i = 0;
    for
    (
        DLList<FoamX::ITypeDescriptorImpl*>::iterator iter =
            subTypes_.begin();
        iter != subTypes_.end();
        ++iter
    )
    {
        // Add reference to the sub type.
        (*typeList)[i++] = iter()->_this();
    }

    // Return the new sub-type list.
    return typeList;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::ITypeDescriptor_ptr FoamX::ITypeDescriptorImpl::elementType()
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::subTypes()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Return a reference to the TypeDescriptor object.
    ITypeDescriptor_ptr pTypeDesc = FoamXServer::ITypeDescriptor::_nil();

    try
    {
        if
        (
            type_ != FoamXServer::Type_List
         && type_ != FoamXServer::Type_FixedList
        )
        {
            throw FoamXError
            (
                E_FAIL,
                "Attempt to call elementType() for non-List type "
              + FoamXTypes::typeName(type_) + " " + name_,
                functionName,
                __FILE__, __LINE__
            );
        }

        if (subTypes_.size() != 1)
        {
            throw FoamXError
            (
                E_FAIL,
                "Number of subTypes != 1 for List type descriptor " + name_,
                functionName,
                __FILE__, __LINE__
            );
        }

        pTypeDesc = FoamXServer::ITypeDescriptor::_duplicate
        (
            subTypes_.first()->_this()
        );
    }
    CATCH_ALL(functionName);

    return pTypeDesc;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ITypeDescriptorImpl::addSubType
(
    FoamXType type,
    ITypeDescriptor_out subEntry
)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::addSubType"
        "(FoamXType type, ITypeDescriptor_out subEntry)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Make sure this is a compound type.
        if (!FoamXTypes::isCompound(type_))
        {
            throw FoamXError
            (
                E_FAIL,
                "Can't add a sub type to a non-compound type "
               + FoamXTypes::typeName(type_) + " " + name_,
                functionName,
                __FILE__, __LINE__
            );
        }

        // Check whether we can add a sub type -
        // Only dictionaries and compound types can have more than one sub-type.
        if
        (
            (type_ == Type_List || type_ == Type_FixedList)
         && subTypes_.size() >= 1
        )
        {
            throw FoamXError
            (
                E_FAIL,
                "One sub-type only for "
                + FoamXTypes::typeName(type_) + " " + name_,
                functionName,
                __FILE__, __LINE__
            );
        }

        // Temporary name for this sub-type.
        word name = "newSubType";
        name += name(subTypes_.size());

        // Create type descriptor object.
        ITypeDescriptorImpl* pTypeDescriptor = new ITypeDescriptorImpl
        (
            name,
            type,
            path_
        );

        // Add to map and the end of the sub type list.
        subTypes_.append(pTypeDescriptor);

        log << "Added sub-type " << name << "." << endl;

        // Return a reference to the new TypeDescriptor object.
        subEntry = pTypeDescriptor->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ITypeDescriptorImpl::removeSubType(ITypeDescriptor_ptr subEntry)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::removeSubType"
        "(ITypeDescriptor_ptr subEntry)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Make sure this is a compound type.
        if (!FoamXTypes::isCompound(type_))
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Can't remove a sub type from a non-compound type "
              + FoamXTypes::typeName(type_) + " " + name_,
                functionName,
                __FILE__, __LINE__
            );
        }

        // Remove from List.
        for
        (
            DLList<ITypeDescriptorImpl*>::iterator iter = 
                subTypes_.begin();
            iter != subTypes_.end();
            ++iter
        )
        {
            ITypeDescriptor_var pTypeDesc = iter()->_this();

            // Check for equivalence.
            if (pTypeDesc->_is_equivalent(subEntry))
            {
                // Release and remove from linked list.
                subTypes_.remove(iter);
                break;
            }
        }

        log << "Removed sub-type." << endl;
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

CORBA::Boolean FoamX::ITypeDescriptorImpl::hasDefaultValue()
{
    return default_ != NULL;
}


void FoamX::ITypeDescriptorImpl::getDefaultValue
(
    IDictionaryEntry_out defaultValue
)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::getDefaultValue"
        "(IDictionaryEntry_out default)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        if (default_ == NULL)
        {
            /*
            // Don't throw error until the Java side checks hasDefaultValue()
            // and requests the creation of the default value only if required
            throw FoamXError
            (
                E_UNEXPECTED,
                "Invalid default object for "
              + FoamXTypes::typeName(type_) + " " + name_,
                functionName,
                __FILE__, __LINE__
            );
            */

            default_ = new IDictionaryEntryImpl(_this());
        }

        defaultValue = default_->_this();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ITypeDescriptorImpl::validate()
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::validate()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool FoamX::ITypeDescriptorImpl::setType(const word& type)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::setType(const word& type)";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Reset all types.
    type_ = Type_Undefined;
    minValue_.setType(Type_Undefined);
    maxValue_.setType(Type_Undefined);

    if (FoamXTypes::found(type))
    {
        type_ = FoamXTypes::lookupType(type);

        // Set the type of the Any parameters.
        if (!isCompoundType())
        {
            if (FoamXTypes::isNumber(type_))
            {
                minValue_.setType(type_);
                maxValue_.setType(type_);
            }

            // Set the type of the values on the value list.
            forAll(valueList_, i)
            {
                valueList_[i].setType(type_);
            }
        }

        return true;
    }
    else
    {
        return false;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ITypeDescriptorImpl::load
(
    const dictionary& typeDict,
    const dictionary& foamTypesDict
)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::load"
        "(const dictionary& typeDict, "
        "const dictionary& foamTypesDict)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        log << "Reading " << typeDict.name()
            << " start line " << typeDict.startLineNumber()
            << " end line " << typeDict.endLineNumber() << endl;

        // Check that the dictionary has the minimum required information.
        if (!typeDict.found("type"))
        {
            throw FoamXError
            (
                E_FAIL,
                "Mandatory entry 'type' not found in dictionary '"
              + typeDict.name() + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Get the type

        // If it is a dictionary get the type from the dictionary
        if (typeDict.isDict("type"))
        {
            // Get the type dictionary ...
            const dictionary& typeTypeDict(typeDict.subDict("type"));

            // ... and check it only has one entry
            if (typeTypeDict.size() != 1)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Dictionary for 'type' of type '"  
                  + FoamXTypes::typeName(type_) + " " + name_
                  + " in dictionary '" + typeDict.name()
                  + "' does contain a single entry",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Get the type entry ...
            const entry& typeEntry = *typeTypeDict.first();

            // ... and it's dictionary ...
            const dictionary& typeEntryTypeDict = typeEntry.dict();

            // ... and get the actual type from it ...
            word type(typeEntryTypeDict.lookup("type"));

            // ... then resolve this type as normal ...
            resolveType(type, typeEntryTypeDict, foamTypesDict);

            // ... but also add the optional additiona entries and data
            addEntries(typeDict, foamTypesDict);
            readOptionalData(typeDict);
        }
        else // it's a system type or defined in foamTypesDict
        {
            // Get the type
            word type(typeDict.lookup("type"));

            // Resolve any Foam types and get sub types for compound types.
            resolveType(type, typeDict, foamTypesDict);
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ITypeDescriptorImpl::resolveType
(
    word& type,
    const dictionary& typeDict,
    const dictionary& foamTypesDict
)
{
    if (typeDict.found("types"))
    {
        resolveTypeDict
        (
            type,
            typeDict,
            foamTypesDict + typeDict.subDict("types")
        );
    }
    else
    {
        resolveTypeDict(type, typeDict, foamTypesDict);
    }
}


void FoamX::ITypeDescriptorImpl::resolveTypeDict
(
    word& type,
    const dictionary& typeDict,
    const dictionary& foamTypesDict
)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::resolveTypeDict(word& type, "
        "const dictionary& typeDict, "
        "const dictionary& foamTypesDict)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        log << "Reading " << typeDict.name()
            << " start line " << typeDict.startLineNumber()
            << " end line " << typeDict.endLineNumber() << endl;

        // See if the specified type is a FoamX type.
        if (FoamXTypes::found(type))
        {
            readType(type, typeDict, foamTypesDict);
        }
        else
        {
            // Specified type may be in the Foam types dictionary.
            if (foamTypesDict.found(type))
            {
                const dictionary& typeAliasDict(foamTypesDict.subDict(type));
                typeAliasDict.lookup("type") >> type;

                resolveType(type, typeAliasDict, foamTypesDict);
                addEntries(typeDict, foamTypesDict);
                readOptionalData(typeDict);
            }
            else
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Type '" + type + "' not defined for "
                  + FoamXTypes::typeName(type_) + " " + name_,
                    functionName,
                    __FILE__, __LINE__
                );
            }
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ITypeDescriptorImpl::readType
(
    word& type,
    const dictionary& typeDict,
    const dictionary& foamTypesDict
)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::readType(word& type, "
        "const dictionary& typeDict, "
        "const dictionary& foamTypesDict)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        log << "Reading " << typeDict.name()
            << " start line " << typeDict.startLineNumber()
            << " end line " << typeDict.endLineNumber() << endl;

        // Set the type. Should have been resolved to a primitive type by now.
        if (!setType(type))
        {
            throw FoamXError
            (
                E_FAIL,
                "Illegal type token '" + type + "' of type '"
              + FoamXTypes::typeName(type_) + " " + name_
              + "' in dictionary '" + typeDict.name() + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // If this type is a compound type, we need the sub type list.
        if (FoamXTypes::isCompound(type_))
        {
            if (type_ == Type_FixedList)
            {
                addElementType(typeDict, foamTypesDict);

                if (!typeDict.found("numElements"))
                {
                    throw FoamXError
                    (
                        E_FAIL,
                        "Mandatory entry 'numElements' of type '" 
                      + FoamXTypes::typeName(type_) + " " + name_
                      + "' not found in dictionary '" + typeDict.name() + "'.",
                        functionName,
                        __FILE__, __LINE__
                    );
                }
                typeDict.lookup("numElements") >> numElements_;

                if (!typeDict.found("elementLabels"))
                {
                    throw FoamXError
                    (
                        E_FAIL,
                        "Mandatory entry 'elementLabels' of type '"  
                      + FoamXTypes::typeName(type_) + " " + name_
                      + "' not found in dictionary '" + typeDict.name() + "'.",
                        functionName,
                        __FILE__, __LINE__
                    );
                }
                elementLabels_.read(typeDict.lookup("elementLabels"));
            }
            else if (type_ == Type_List)
            {
                addElementType(typeDict, foamTypesDict);
            }
            else
            {
                addEntries(typeDict, foamTypesDict);
            }
        }

        readOptionalData(typeDict);
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ITypeDescriptorImpl::addEntries
(
    const dictionary& typeDict,
    const dictionary& foamTypesDict
)
{
    if
    (
        (
            type_ == Type_Dictionary
         || type_ == Type_Selection
         || type_ == Type_Compound
        )
     && typeDict.found("entries")
    )
    {
        if (type_ == Type_Dictionary)
        {
            addDictionaryEntries(typeDict, foamTypesDict);
        }
        else
        {
            addCompoundEntries(typeDict, foamTypesDict);
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ITypeDescriptorImpl::addDictionaryEntries
(
    const dictionary& typeDict,
    const dictionary& foamTypesDict
)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::addDictionaryEntries"
        "(const dictionary& typeDict, const dictionary& foamTypesDict)";

    const dictionary& entries(typeDict.subDict("entries"));

    for
    (
        dictionary::const_iterator iter = entries.begin();
        iter != entries.end();
        ++iter
    )
    {
        const word& keyword = iter().keyword();

        if (iter().isDict())
        {
            subTypes_.append
            (
                new ITypeDescriptorImpl
                (
                    keyword,
                    path_,
                    iter().dict(),
                    foamTypesDict
                )
            );
        }
        else
        {
            if (!iter().stream().size())
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Entry '" + keyword
                  + "' for 'entries' of type '"  
                  + FoamXTypes::typeName(type_) + " " + name_
                  + "' in dictionary '" + typeDict.name()
                  + "' does not have a type.",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            word subType(iter().stream());

            if (foamTypesDict.found(subType))
            {
                subTypes_.append
                (
                    new ITypeDescriptorImpl
                    (
                        keyword,
                        path_,
                        foamTypesDict.subDict(subType),
                        foamTypesDict
                    )
                );
            }
            else if (FoamXTypes::found(subType))
            {
                subTypes_.append
                (
                    new ITypeDescriptorImpl
                    (
                        keyword,
                        FoamXTypes::lookupType(subType),
                        path_
                    )
                );
            }
            else
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Type '" + subType
                  + "' for 'entries' of type '"  
                  + FoamXTypes::typeName(type_) + " " + name_
                  + "' in dictionary '" + typeDict.name() + "' not defined.",
                    functionName,
                    __FILE__, __LINE__
                );
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ITypeDescriptorImpl::addCompoundEntries
(
    const dictionary& typeDict,
    const dictionary& foamTypesDict
)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::addCompoundEntries"
        "(const dictionary& typeDict, const dictionary& foamTypesDict)";

    const dictionary& entries(typeDict.subDict("entries"));

    for
    (
        dictionary::const_iterator iter = entries.begin();
        iter != entries.end();
        ++iter
    )
    {
        const word& subType = iter().keyword();

        if (iter().isDict())
        {
            subTypes_.append
            (
                new ITypeDescriptorImpl
                (
                    subType,
                    path_,
                    iter().dict(),
                    foamTypesDict
                )
            );
        }
        else if (foamTypesDict.found(subType))
        {
            subTypes_.append
            (
                new ITypeDescriptorImpl
                (
                    subType,
                    path_,
                    foamTypesDict.subDict(subType),
                    foamTypesDict
                )
            );
        }
        else if (FoamXTypes::found(subType))
        {
            subTypes_.append
            (
                new ITypeDescriptorImpl
                (
                    subType,
                    FoamXTypes::lookupType(subType),
                    path_
                )
            );
        }
        else
        {
            throw FoamXError
            (
                E_FAIL,
                "Type '" + subType
              + "' for 'entries' of type '"  
              + FoamXTypes::typeName(type_) + " " + name_
              + "' in dictionary '" + typeDict.name() + "' not defined.",
                functionName,
                __FILE__, __LINE__
            );
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ITypeDescriptorImpl::addElementType
(
    const dictionary& typeDict,
    const dictionary& foamTypesDict
)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::addSubType(const dictionary& typeDict, "
        "const dictionary& foamTypesDict)";

    LogEntry log(functionName, __FILE__, __LINE__);

    if (!typeDict.found("elementType"))
    {
        throw FoamXError
        (
            E_FAIL,
            "Mandatory entry 'elementType' of type '"  
          + FoamXTypes::typeName(type_) + "' '" + name_
          + "' not found in dictionary '" + typeDict.name() + "'.",
            functionName,
            __FILE__, __LINE__
        );
    }

    if (typeDict.isDict("elementType"))
    {
        const dictionary& elementTypeDict = typeDict.subDict("elementType");

        if (elementTypeDict.size() != 1)
        {
            throw FoamXError
            (
                E_FAIL,
                "Dictionary for 'elementType' of type '"  
              + FoamXTypes::typeName(type_) + " " + name_
              + " in dictionary '" + typeDict.name()
              + "' does contain a single entry",
                functionName,
                __FILE__, __LINE__
            );
        }

        const entry& elementTypeEntry = *elementTypeDict.first();

        subTypes_.append
        (
            new ITypeDescriptorImpl
            (
                elementTypeEntry.keyword(),
                path_,
                elementTypeEntry.dict(),
                foamTypesDict
            )
        );
    }
    else
    {
        word subType = typeDict.lookup("elementType");

        if (FoamXTypes::found(subType))
        {
            // Add to map and list.
            subTypes_.append
            (
                new ITypeDescriptorImpl
                (
                    subType,
                    FoamXTypes::lookupType(subType),
                    path_
                )
            );
        }
        else if (foamTypesDict.found(subType))
        {
            subTypes_.append
            (
                new ITypeDescriptorImpl
                (
                    subType,
                    path_,
                    foamTypesDict.subDict(subType),
                    foamTypesDict
                )
            );
        }
        else
        {
            throw FoamXError
            (
                E_FAIL,
                "Type '" + subType
              + "' for 'elementType' of type '"  
              + FoamXTypes::typeName(type_) + " " + name_
              + " in dictionary '" + typeDict.name() + "' not defined.",
                functionName,
                __FILE__, __LINE__
            );
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ITypeDescriptorImpl::readOptionalData
(
    const dictionary& typeDict
)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::readOptionalData"
        "(const dictionary& typeDict)";

    LogEntry log(functionName, __FILE__, __LINE__);

    log << "Reading " << typeDict.name()
        << " start line " << typeDict.startLineNumber()
        << " end line " << typeDict.endLineNumber() << endl;

    // Read the rest of the common type information.
    if (typeDict.found("displayName"))
    {
        typeDict.lookup("displayName") >> displayName_;
    }

    if (typeDict.found("description"))
    {
        typeDict.lookup("description") >> description_;
    }

    if (typeDict.found("dictionaryPath"))
    {
        typeDict.lookup("dictionaryPath") >> dictionaryPath_;
        dictionaryPath_.expand();
    }

    if (typeDict.found("category"))
    {
        typeDict.lookup("category")  >> category_;
    }

    if (typeDict.found("comment"))
    {
        typeDict.lookup("comment") >> comment_;
    }

    if (typeDict.found("helpURL"))
    {
        typeDict.lookup("helpURL") >> helpURL_;
    }

    if (typeDict.found("iconURL"))
    {
        typeDict.lookup("iconURL") >> iconURL_;
    }

    if (typeDict.found("optional"))
    {
        optional_ = readBool(typeDict.lookup("optional"));
    }

    if (typeDict.found("visible"))
    {
        visible_ = readBool(typeDict.lookup("visible"));
    }

    if (typeDict.found("editable"))
    {
        editable_ = readBool(typeDict.lookup("editable"));
    }


    // Read non-compound specific information.
    if (FoamXTypes::isPrimitive(type_))
    {
        if (FoamXTypes::isNumber(type_))
        {
            if (typeDict.found("minValue"))
            {
                minValue_.read(typeDict.lookup("minValue"));
            }
            else
            {
                minValue_.setType(type_);
                minValue_.setMin();
            }

            if (typeDict.found("maxValue"))
            {
                maxValue_.read(typeDict.lookup("maxValue"));
            }
            else
            {
                maxValue_.setType(type_);
                maxValue_.setMax();
            }
        }

        if (typeDict.found("lookupDict"))
        {
            typeDict.lookup("lookupDict") >> lookupDict_;
        }

        // Read value list.
        if (typeDict.found("valueList"))
        {
            Istream& is = typeDict.lookup("valueList");

            is.readBeginList("List");

            SLList<FoamXAny> sllist;

            // Loop over list elements until ')'
            token lastToken(is);
            while
            (
                !(
                    lastToken.isPunctuation()
                 && lastToken.pToken() == token::END_LIST
                )
            )
            {
                is.putBack(lastToken);
                sllist.append(FoamXAny(type_, is));
                is >> lastToken;
            }

            valueList_ = sllist;
        }
    }


    if (typeDict.found("options"))
    {
        const dictionary& subTypeDefaults(typeDict.subDict("options"));

        for
        (
            DLList<FoamX::ITypeDescriptorImpl*>::iterator iter =
                subTypes_.begin();
            iter != subTypes_.end();
            ++iter
        )
        {
            if (subTypeDefaults.found(iter()->name()))
            {
                iter()->readOptionalData
                (
                    subTypeDefaults.subDict(iter()->name())
                );
            }
        }
    }

    if (typeDict.found("default"))
    {
        if (!default_)
        {
            default_ = new IDictionaryEntryImpl(_this());
        }

        default_->load(typeDict.lookupEntry("default"));

        DLList<FoamX::ITypeDescriptorImpl*>::iterator iter1 =
            subTypes_.begin();

        DLList<IDictionaryEntryImpl*>::const_iterator iter2 =
            default_->subElementPtrs().begin();

        for
        (
            ;
            iter1 != subTypes_.end()
         && iter2 != default_->subElementPtrs().end();
            ++iter1, ++iter2
        )
        {
            if (!iter1()->default_)
            {
                iter1()->default_ = new IDictionaryEntryImpl(iter1()->_this());
            }

            *iter1()->default_ = *iter2();
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ITypeDescriptorImpl::save
(
    FoamX::DictionaryWriter& dictWriter,
    bool topLevelDict
)
{
    static const char* functionName =
        "FoamX::ITypeDescriptorImpl::save"
        "(FoamX::DictionaryWriter& dictWriter, bool topLevelDict)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        if (!topLevelDict)
        {
            // Start a sub-dictionary for this type.
            dictWriter.startSubDict(name_);
        }
        else
        {
            dictWriter.writeEndl();
        }

        // Write type.
        dictWriter.writeEntry("type", FoamXTypes::typeName(type_));

        // Write optional common data.
        if (displayName_.length() > 0)
        {
            dictWriter.writeEntry("displayName", displayName_);
        }

        if (description_.length() > 0)
        {
            dictWriter.writeEntry("description", description_);
        }
        
        if (category_.length() > 0)
        {
            dictWriter.writeEntry("category", category_);
        }

        if (comment_.length() > 0)
        {
            dictWriter.writeEntry("comment", comment_);
        }

        if (helpURL_.length() > 0)
        {
            dictWriter.writeEntry("helpURL", helpURL_);
        }

        if (iconURL_.length() > 0)
        {
            dictWriter.writeEntry("iconURL", iconURL_);
        }

        // Only write optional flag if true.
        if (optional_)
        {
            dictWriter.writeEntry("optional", optional_);
        }

        if (visible_)
        {
            dictWriter.writeEntry("visible", visible_);
        }

        if (editable_)
        {
            dictWriter.writeEntry("editable", editable_);
        }

        // Write non-compound data.
        if (FoamXTypes::isPrimitive(type_))
        {
            if (FoamXTypes::isNumber(type_))
            {
                if (minValue_ != FoamXAny(type_))
                {
                    dictWriter.writeEntry("minValue", minValue_);
                }

                if (maxValue_ != FoamXAny(type_))
                {
                    dictWriter.writeEntry("maxValue", maxValue_);
                }
            }

            if (lookupDict_.length() > 0)
            {
                dictWriter.writeEntry("lookupDict", lookupDict_);
            }

            if (valueList_.size() > 0)
            {
                dictWriter.writeEntry("valueList", valueList_);
            }
        }
        else if(FoamXTypes::isCompound(type_))
        {
            if (type_ == Type_FixedList)
            {
                dictWriter.writeEntry("numElements", numElements_);

                dictWriter.writeEntry
                (
                    "elementType",
                    word(subTypes_.first()->name())
                );

                dictWriter.writeEntry("elementLabels", elementLabels_);
            }
            else if (type_ == Type_List)
            {
                dictWriter.writeEntry
                (
                    "elementType",
                    word(subTypes_.first()->name())
                );
            }
            else
            {
                // Write compound data.
                if (type_ == Type_Dictionary)
                {
                    if (dictionaryPath_.length() > 0)
                    {
                        dictWriter.writeEntry
                        (
                            "dictionaryPath",
                            dictionaryPath_
                        );
                    }
                }

                dictWriter.writeKeyword("entries");
                dictWriter.startList(subTypes_.size());

                label i = 0;
                // Write sub-type data into their sub-dictionaries.
                for
                (
                    DLList<ITypeDescriptorImpl*>::iterator iter = 
                        subTypes_.begin();
                    iter != subTypes_.end();
                    ++iter
                )
                {
                    iter()->save(dictWriter, false);

                    if (++i < subTypes_.size())
                    {
                        dictWriter.writeEndl();
                    }
                }

                dictWriter.endList();
                dictWriter.endEntry();
            }
        }

        // Validate the default_ dictionary entry items.
        if (default_)
        {
            default_->validate();

            dictWriter.writeKeyword("default");
            default_->save(dictWriter, false);
            dictWriter.endEntry();
        }

        // Close this sub-dictionary if required.
        if (!topLevelDict)
        {
            dictWriter.endSubDict();
        }
    }
    CATCH_ALL(functionName);
}


// ************************************************************************* //
