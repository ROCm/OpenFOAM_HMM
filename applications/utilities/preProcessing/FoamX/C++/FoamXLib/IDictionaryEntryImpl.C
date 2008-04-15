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
#include "vector.H"
#include "tensor.H"
#include "typeInfo.H"

// Project header files.
#include "IDictionaryEntryImpl.H"
#include "FoamX.H"
#include "FoamXErrors.H"
#include "LogEntry.H"
#include "ITypeDescriptorImpl.H"

// Namespaces
#include "FoamXNameSpaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IDictionaryEntryImpl::IDictionaryEntryImpl
(
    ITypeDescriptor_ptr typeDesc
)
:
    modified_(false),
    listTokenPtr_(NULL),
    selection_(0)
{
    static const char* functionName =
        "FoamX::IDictionaryEntryImpl::IDictionaryEntryImpl"
        "(ITypeDescriptor_ptr typeDesc)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Check parameters.
        if (CORBA::is_nil(typeDesc))
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Invalid TypeDescriptor reference",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Bind the given field type descriptor.
        bindType(typeDesc);
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IDictionaryEntryImpl::IDictionaryEntryImpl
(
    IDictionaryEntryImpl& dictEntry
)
:
    PortableServer::ServantBase(),
    PortableServer::StaticImplementation(),
    POA_FoamXServer::IDictionaryEntry(),
    PortableServer::RefCountServantBase(),
    modified_(false),
    listTokenPtr_(NULL),
    selection_(0)
{
    static const char* functionName =
        "FoamX::IDictionaryEntryImpl::IDictionaryEntryImpl"
        "(ITypeDescriptor_ptr typeDesc)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Bind the given field type descriptor.
        bindType(dictEntry.typeDescriptor());

        operator=(dictEntry);
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IDictionaryEntryImpl::~IDictionaryEntryImpl()
{
    static const char* functionName =
        "FoamX::IDictionaryEntryImpl::~IDictionaryEntryImpl()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Release all sub-element objects.
        for
        (
            DLList<IDictionaryEntryImpl*>::iterator iter =
                subElements_.begin();
            iter != subElements_.end();
            ++iter
        )
        {
            iter()->_remove_ref();
        }
        subElements_.clear();

        if (listTokenPtr_)
        {
            delete listTokenPtr_;
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::ITypeDescriptor_ptr FoamX::IDictionaryEntryImpl::typeDescriptor()
{
    static const char* functionName =
        "FoamX::IDictionaryEntryImpl::typeDescriptor()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Return a reference to the TypeDescriptor object.
    ITypeDescriptor_ptr pTypeDesc = ITypeDescriptor::_nil();

    if (!CORBA::is_nil(typeDescriptor_))
    {
        pTypeDesc = ITypeDescriptor::_duplicate(typeDescriptor_);
    }

    return pTypeDesc;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::FoamXAny* FoamX::IDictionaryEntryImpl::value()
{
    static const char* functionName =
        "FoamX::IDictionaryEntryImpl::value()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Duplicate and return.
    return new FoamXServer::FoamXAny(value_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IDictionaryEntryImpl::value(const FoamXServer::FoamXAny& newValue)
{
    static const char* functionName =
        "FoamX::IDictionaryEntryImpl::value(const CORBA::Any& newValue)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        if (typeDescriptor_->editable() && value_ != newValue)
        {
            // Set the value of the any object.
            // Will throw an exception if the incoming any is not compatible.
            value_.setValue(newValue);

            // Set modified flag.
            modified_ = true;
        }
    }
    catch (...)
    {}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IDictionaryEntryImpl::setValue
(
    const FoamXServer::FoamXAny& newValue
)
{
    static const char* functionName =
        "FoamX::IDictionaryEntryImpl::setValue(const CORBA::Any& newValue)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        if (typeDescriptor_->editable() && value_ != newValue)
        {
            if (FoamXTypes::isNumber(typeDescriptor_->type()))
            {
                //***HGWnewValue > *typeDescriptor_->maxValue();
                //***HGWnewValue < *typeDescriptor_->minValue();
            }

            // Set the value of the any object.
            // Will throw an exception if the incoming any is not compatible.
            value_.setValue(newValue);

            // Set modified flag.
            modified_ = true;
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class CorbaType>
void FoamX::IDictionaryEntryImpl::expandPrimitiveList()
{
    static const char* functionName =
        "FoamX::IDictionaryEntryImpl::expandPrimitiveList<Type, CorbaType>()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        if (!isA<List<Type> >(listTokenPtr_->compoundToken()))
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "List is not of the expected type",
                functionName,
                __FILE__, __LINE__
            );
        }

        const List<Type>& sl = dynamicCast<const List<Type> >
        (
            listTokenPtr_->compoundToken()
        );

        forAll (sl, i)
        {
            IDictionaryEntryImpl* pSubEntry = new IDictionaryEntryImpl
            (
                typeDescriptor_->elementType()
            );
                
            if (pSubEntry == NULL)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Failed to create "
                    "IDictionaryEntryImpl object",
                    functionName,
                    __FILE__, __LINE__
                );
            }
            
            pSubEntry->value_.value <<= CorbaType(sl[i]);
                
            subElements_.append(pSubEntry);
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Form>
void FoamX::IDictionaryEntryImpl::expandFixedListList()
{
    static const char* functionName =
        "FoamX::IDictionaryEntryImpl::"
        "expandFixedListList<Form, Cmpt, nCmpt>()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        if (!isA<List<Form> >(listTokenPtr_->compoundToken()))
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "List is not of the expected type",
                functionName,
                __FILE__, __LINE__
            );
        }

        const List<Form>& vl = dynamicCast<const List<Form> >
        (
            listTokenPtr_->compoundToken()
        );

        forAll (vl, i)
        {
            IDictionaryEntryImpl* pSubEntry = new IDictionaryEntryImpl
            (
                typeDescriptor_->elementType()
            );
                
            if (pSubEntry == NULL)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Failed to create "
                    "IDictionaryEntryImpl object",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            if (pSubEntry->subElements_.size() != Form::nComponents)
            {
                throw FoamXError
                (
                    E_UNEXPECTED,
                    "Number of sub-elements in Type_FixedList = "
                  + name(pSubEntry->subElements_.size()) + " is not equal to "
                  + "the number of elements in the given Form = "
                  + name(Form::nComponents),
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Loop over all vector-space elements.
            direction j=0;
            for
            (
                DLList<IDictionaryEntryImpl*>::iterator
                    iter = pSubEntry->subElements_.begin();
                iter != pSubEntry->subElements_.end();
                ++iter
            )
            {
                iter()->value_.value <<= vl[i][j++];
            }
                
            subElements_.append(pSubEntry);
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IDictionaryEntryImpl::expandList()
{
    static const char* functionName =
        "FoamX::IDictionaryEntryImpl::expandList()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        if (typeDescriptor_->type() != Type_List)
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Unexpected call to expandList for non-list type "
              + word(typeDescriptor_->name()),
                functionName,
                __FILE__, __LINE__
            );
        }

        if (!listTokenPtr_)
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Unexpected call to expandList for expanded list",
                functionName,
                __FILE__, __LINE__
            );
        }

        if (listTokenPtr_->compoundToken().size() > 1000)
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Unexpected call to expandList for list "
                "larger than 1000 elements",
                functionName,
                __FILE__, __LINE__
            );
        }

        if
        (
            isA<List<label> >(listTokenPtr_->compoundToken())
        )
        {
            expandPrimitiveList<label, CORBA::Long>();
        }
        else if
        (
            isA<List<scalar> >(listTokenPtr_->compoundToken())
        )
        {
            expandPrimitiveList<scalar, CORBA::Double>();
        }
        else if
        (
            isA<List<vector> >(listTokenPtr_->compoundToken())
        )
        {
            expandFixedListList<vector>();
        }
        else if
        (
            isA<List<tensor> >(listTokenPtr_->compoundToken())
        )
        {
            expandFixedListList<tensor>();
        }

        delete listTokenPtr_;
        listTokenPtr_ = NULL;
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::DictionaryEntryList* FoamX::IDictionaryEntryImpl::subElements()
{
    static const char* functionName =
        "FoamX::IDictionaryEntryImpl::subElements()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        if (listTokenPtr_)
        {
            expandList();
        }

        // Scan for all the visible entries
        int nVisibleElements = 0;
        for
        (
            DLList<IDictionaryEntryImpl*>::iterator iter =
                subElements_.begin();
            iter != subElements_.end();
            ++iter
        )
        {
            if (iter()->typeDescriptor_->visible())
            {
                nVisibleElements++;
            }
        }

        // Construct a new DictionaryEntryList.
        DictionaryEntryList* entryList = new DictionaryEntryList();
        entryList->length(nVisibleElements);

        int i = 0;
        for
        (
            DLList<IDictionaryEntryImpl*>::iterator iter =
                subElements_.begin();
            iter != subElements_.end();
            ++iter
        )
        {
            if (iter()->typeDescriptor_->visible())
            {
                // Add a reference to the sub element.
                (*entryList)[i++] = iter()->_this();
            }
        }

        // Return the entry list.
        return entryList;
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

CORBA::Long FoamX::IDictionaryEntryImpl::nSubElements()
{
    if (listTokenPtr_)
    {
        return listTokenPtr_->compoundToken().size();
    }
    else
    {
        return subElements_.size();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

CORBA::Boolean FoamX::IDictionaryEntryImpl::packedList()
{
    if (listTokenPtr_)
    {
        return true;
    }
    else
    {
        return false;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//FoamX::DictionaryEntry* FoamX::IDictionaryEntryImpl::lookup
//(
//    const char* entryName
//) const
//{
//    static const char* functionName =
//        "FoamX::IDictionaryEntryImpl::lookup(const char*)";
//
//    LogEntry log(functionName, __FILE__, __LINE__);
//
//    try
//    {
// 
//        // Construct a new DictionaryEntryList.
//        DictionaryEntry* entry = new DictionaryEntry();
//
//        int i = 0;
//        for
//        (
//            DLList<IDictionaryEntryImpl*>::iterator iter =
//                subElements_.begin();
//            iter != subElements_.end();
//            ++iter
//        )
//        {
//            if (iter()->typeDescriptor_->name() == entryName)
//            {
//                return new DictionaryEntry(iter()->_this());
//            }
//        }
//
//        throw FoamXError
//        (
//            E_FAIL,
//            "Did not find " + entryName,
//            functionName,
//            __FILE__, __LINE__
//        );
//    }
//    CATCH_ALL(functionName);
//}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

CORBA::Long FoamX::IDictionaryEntryImpl::selection()
{
    return selection_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IDictionaryEntryImpl::selection(CORBA::Long newSelection)
{
    if (selection_ != newSelection)
    {
        selection_ = newSelection;

        // Set modified flag.
        modified_ = true;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IDictionaryEntryImpl::addElement(IDictionaryEntry_out subEntry)
{
    static const char* functionName =
        "FoamX::IDictionaryEntryImpl::addElement"
        "(IDictionaryEntry_out subEntry)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        if (typeDescriptor_->editable())
        {
            // Make sure this is a list.
            if (typeDescriptor_->type() != Type_List)
            {
                throw FoamXError
                (
                    E_UNEXPECTED,
                    "Unexpected call to addElement for non-list type "
                  + word(typeDescriptor_->name()),
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Create an appropriate sub entry object and store its reference.
            IDictionaryEntryImpl* pSubEntry = new IDictionaryEntryImpl
            (
                typeDescriptor_->elementType()
            );

            if (pSubEntry == NULL)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Couldn't create IDictionaryEntryImpl object for "
                  + word(typeDescriptor_->name()),
                    functionName,
                    __FILE__, __LINE__
                );
            }
            subElements_.append(pSubEntry);

            // Return a reference to the new element.
            subEntry = pSubEntry->_this();

            // Set modified flag.
            modified_ = true;
        }
        else
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Unexpected call to addElement for non-editable type "
              + word(typeDescriptor_->name()),
                functionName,
                __FILE__, __LINE__
            );
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IDictionaryEntryImpl::removeElement
(
    IDictionaryEntry_ptr subEntry
)
{
    static const char* functionName =
        "FoamX::IDictionaryEntryImpl::removeElement"
        "(IDictionaryEntry_ptr subEntry)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        if (typeDescriptor_->editable())
        {
            // Make sure this is a list.
            if (typeDescriptor_->type() != Type_List)
            {
                throw FoamXError
                (
                    E_UNEXPECTED,
                    "Unexpected call to removeElement for non-list type "
                  + word(typeDescriptor_->name()),
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Remove from List.
            for
            (
                DLList<IDictionaryEntryImpl*>::iterator iter = 
                    subElements_.begin();
                iter != subElements_.end();
                ++iter
            )
            {
                IDictionaryEntry_var pElement = iter()->_this();

                // Check for equivalence.
                if (pElement->_is_equivalent(subEntry))
                {
                    // Release and remove from linked list.
                    iter()->_remove_ref();
                    subElements_.remove(iter);
                    break;
                }
            }

            // Set modified flag.
            modified_ = true;
        }
        else
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Unexpected call to removeElement for non-editable type "
              + word(typeDescriptor_->name()),
                functionName,
                __FILE__, __LINE__
            );
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IDictionaryEntryImpl::clearSubElements()
{
    static const char* functionName =
        "FoamX::IDictionaryEntryImpl::clearSubElements()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        if
        (
            static_cast<ITypeDescriptor*>(typeDescriptor_)
         == ITypeDescriptor::_nil()
         || typeDescriptor_->editable()
        )
        {
            // Release all sub-element objects.
            for
            (
                DLList<IDictionaryEntryImpl*>::iterator iter =
                    subElements_.begin();
                iter != subElements_.end();
                ++iter
            )
            {
                iter()->_remove_ref();
            }
            subElements_.clear();
        }
        else
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Unexpected call to clearSubElement for non-editable type "
              + word(typeDescriptor_->name()),
                functionName,
                __FILE__, __LINE__
            );
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IDictionaryEntryImpl::validate()
{
    static const char* functionName =
        "FoamX::IDictionaryEntryImpl::validate()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // If primitive entry, make sure that we have a valid value
        if (!typeDescriptor_->isPrimitiveType())
        {
            // If this entry is not optional,
            // make sure that a value has been set.
            if (!value_.IsSet() && !typeDescriptor_->optional())
            {
                // Return the path of this entry.
                throw ValidationError
                (
                    E_INVALID_ARG,
                    "Non-optional entry does not have a value",
                    typeDescriptor_->path()
                );
            }

            // Check that the value is in range.
            if (FoamXTypes::isNumber(typeDescriptor_->type()))
            {
                value_ > *typeDescriptor_->maxValue();
                value_ < *typeDescriptor_->minValue();
            }
        }
        else
        {
            // Check all sub-entries.
            for
            (
                DLList<IDictionaryEntryImpl*>::iterator iter =
                    subElements_.begin();
                iter != subElements_.end();
                ++iter
            )
            {
                iter()->validate();
            }
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

CORBA::Boolean FoamX::IDictionaryEntryImpl::modified()
{
    static const char* functionName =
        "FoamX::IDictionaryEntryImpl::modified()";

    LogEntry log(functionName, __FILE__, __LINE__);

    if (modified_) return true;

    try
    {
        // Check all sub-entries.
        for
        (
            DLList<IDictionaryEntryImpl*>::iterator iter =
                subElements_.begin();
            iter != subElements_.end();
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

void FoamX::IDictionaryEntryImpl::save()
{
    // Can't save individual entries.
    // Overridden for root dictionary entries.
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IDictionaryEntryImpl::load
(
    const entry& ent
)
{
    static const char* functionName =
        "FoamX::IDictionaryEntryImpl::load"
        "(const entry&, bool allowNonOptional)";

    if (ent.isStream())
    {
        load(ent.stream());
    }
    else if (ent.isDict())
    {
        load(ent.dict());
    }
    else
    {
        throw FoamXError
        (
            E_UNEXPECTED,
            "Entry is neither a stream nor a dictionary",
            functionName,
            __FILE__, __LINE__
        );
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IDictionaryEntryImpl::load
(
    const dictionary& dict,
    bool allowNonOptional
)
{
    static const char* functionName =
        "FoamX::IDictionaryEntryImpl::load"
        "(const dictionary&, bool allowNonOptional)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        if (typeDescriptor_->type() != Type_Dictionary)
        {
            throw FoamXIOError
            (
                "IDictionaryEntryImpl " + word(typeDescriptor_->name())
              + " is not a dictionary\n",
                dict.name(), dict.startLineNumber(), dict.endLineNumber(),
                functionName,
                __FILE__, __LINE__
            );
        }

        log << "Reading Dictionary" << endl;

        // Loop over all elements.
        for
        (
            DLList<IDictionaryEntryImpl*>::iterator iter
                = subElements_.begin();
            iter != subElements_.end();
            ++iter
        )
        {
            // Get the sub-entry object.
            IDictionaryEntryImpl* pSubElement = iter();
            
            // Get the TypeDescriptor for this entry.
            ITypeDescriptor_var pSubTypeDesc = pSubElement->typeDescriptor_;
                
            // Get the entry name.
            word paramName(pSubTypeDesc->name());
                
            // Make sure this entry exists in the dictionary.
            if (dict.found(paramName))
            {
                pSubElement->load(dict.lookupEntry(paramName));
            }
            else if (!allowNonOptional && !pSubTypeDesc->optional())
            {
                throw FoamXIOError
                (
                    "Non-optional dictionary entry '" + paramName
                  + "' not found in dictionary\n" + dict.name(),
                    dict.name(), dict.startLineNumber(), dict.endLineNumber(),
                    functionName,
                    __FILE__, __LINE__
                );
            }
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IDictionaryEntryImpl::load(Istream& is)
{
    static const char* functionName =
        "FoamX::IDictionaryEntryImpl::load(Istream&)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        log << "Loading DictionaryEntry for "
            << CORBA::String_var(typeDescriptor_->path()) << endl;

        log << "Reading " << is.name() << " line " << is.lineNumber() << endl;

        // See if this is a primitive type.
        if (typeDescriptor_->isPrimitiveType())
        {
            log << "Reading non-compound value" << endl;

            // This entry is not compound so read the single value.
            // Will throw an exception if a value of the required type
            // cannot be extracted from the Istream.
            value_.read(is);
        }
        else
        {
            switch(typeDescriptor_->type())
            {
                case Type_FixedList:
                {
                    log << "Reading FixedList" << endl;

                    // Read beginning of list contents.
                    is.readBeginList("List");
                
                    // Loop over all vector elements.
                    for
                    (
                        DLList<IDictionaryEntryImpl*>::iterator
                            iter = subElements_.begin();
                        iter != subElements_.end();
                        ++iter
                    )
                    {
                        // Get the sub-entry object.
                        IDictionaryEntryImpl* pSubElement = iter();
                
                        // Load the element values from the Istream.
                        pSubElement->load(is);
                    }
                
                    // Read end of contents.
                    is.readEndList("List");
                }
                break;

                case Type_List:
                {
                    log << "Reading List" << endl;
                
                    // Clear old list elements since size of
                    // list may have changed.
                    clearSubElements();
                
                    token firstToken(is);

                    is.fatalCheck
                    (
                        "IDictionaryEntryImpl::load(Istream&, bool)"
                        " : reading first token"
                    );

                    if (firstToken.isCompound())
                    {
                        if (listTokenPtr_)
                        {
                            delete listTokenPtr_;
                        }
                        listTokenPtr_ = new token(firstToken);
                    }
                    else if (firstToken.isLabel())
                    {
                        label s = firstToken.labelToken();

                        // Read beginning of contents
                        char listDelimiter = is.readBeginList("Type_List");

                        if (s)
                        {
                            if (listDelimiter == token::BEGIN_LIST)
                            {
                                for (register label i=0; i<s; i++)
                                {
                                    // Create an appropriate sub entry object
                                    // and store its reference.
                                    IDictionaryEntryImpl* pSubEntry = 
                                    new IDictionaryEntryImpl
                                    (
                                        typeDescriptor_->elementType()
                                    );

                                    if (pSubEntry == NULL)
                                    {
                                        throw FoamXError
                                        (
                                            E_FAIL,
                                            "Failed to create "
                                            "IDictionaryEntryImpl object",
                                            functionName,
                                            __FILE__, __LINE__
                                        );
                                    }
                
                                    // Load the values from the Istream.
                                    pSubEntry->load(is);
                
                                    // Add to list.
                                    subElements_.append(pSubEntry);
                                }
                            }
                            else
                            {
                                // Create an appropriate sub entry object
                                // and store its reference.
                                IDictionaryEntryImpl* pSubEntry =
                                new IDictionaryEntryImpl
                                (
                                    typeDescriptor_->elementType()
                                );

                                if (pSubEntry == NULL)
                                {
                                    throw FoamXError
                                    (
                                        E_FAIL,
                                        "Failed to create IDictionaryEntryImpl "
                                        "object",
                                        functionName,
                                        __FILE__, __LINE__
                                    );
                                }
                
                                // Load the values from the Istream.
                                pSubEntry->load(is);
                
                                // Add to list.
                                subElements_.append(pSubEntry);
                                
                                for (register label i=1; i<s; i++)
                                {
                                    subElements_.append
                                    (
                                        new IDictionaryEntryImpl(*pSubEntry)
                                    );
                                }
                            }
                        }
                        
                        // Read end of contents
                        is.readEndList("Type_List");
                    }
                    else if (firstToken.isPunctuation())
                    {
                        if (firstToken.pToken() != token::BEGIN_LIST)
                        {
                            FatalIOErrorIn
                            (
                                "IDictionaryEntryImpl::load"
                                "(Istream&, bool)",
                                is
                            )   << "incorrect first token, '(', found " 
                                << firstToken.info()
                                << exit(FatalIOError);
                        }

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

                            // Create an appropriate sub entry object
                            // and store its reference.
                            IDictionaryEntryImpl* pSubEntry = new 
                            IDictionaryEntryImpl
                            (
                                typeDescriptor_->elementType()
                            );

                            if (pSubEntry == NULL)
                            {
                                throw FoamXError
                                (
                                    E_FAIL,
                                    "Failed to create "
                                    "IDictionaryEntryImpl object",
                                    functionName,
                                    __FILE__, __LINE__
                                );
                            }
                
                            // Load the values from the Istream.
                            pSubEntry->load(is);
                
                            // Add to list.
                            subElements_.append(pSubEntry);

                            is >> lastToken;
                        }
                    }
                    else
                    {
                        FatalIOErrorIn
                        (
                            "IDictionaryEntryImpl::load(Istream&, bool)",
                            is
                        )   << "incorrect first token, "
                              "expected <int> or '(', found "
                            << firstToken.info()
                            << exit(FatalIOError);
                    }
                }
                break;

                case Type_Dictionary:
                {
                    load(dictionary(is));
                }
                break;

                case Type_Selection:
                {
                    log << "Reading Selection" << endl;

                    // Read selection keyword
                    word keyword(is);

                    log << "Keyword = " << keyword << endl;

                    label i = 0;

                    // Loop over all types to find selection.
                    for
                    (
                        DLList<IDictionaryEntryImpl*>::iterator iter
                            = subElements_.begin();
                        iter != subElements_.end();
                        ++iter, ++i
                    )
                    {
                        // Get the sub-entry object.
                        IDictionaryEntryImpl* pSubElement = iter();

                        if (keyword == pSubElement->typeDescriptor_->name())
                        {
                            selection_ = i;

                            // Load the element values from the Istream.
                            pSubElement->load(is);
                        }
                    }
                }
                break;

                case Type_Compound:
                {
                    log << "Reading Compound" << endl;
                
                    // Loop over all vector elements.
                    for
                    (
                        DLList<IDictionaryEntryImpl*>::iterator iter
                            = subElements_.begin();
                        iter != subElements_.end();
                        ++iter
                    )
                    {
                        // Get the sub-entry object.
                        IDictionaryEntryImpl* pSubElement = iter();
                
                        // Load the element values from the Istream.
                        pSubElement->load(is);
                    }
                }
                break;

                default:
                    throw FoamXError
                    (
                        E_INVALID_ARG,
                        "Invalid Type",
                        functionName,
                        __FILE__, __LINE__
                    );
            }
        }

        // Check state of Istream.
        is.check(functionName);
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IDictionaryEntryImpl::save
(
    FoamX::DictionaryWriter& dictWriter,
    bool dictEntry
)
{
    static const char* functionName =
        "FoamX::IDictionaryEntryImpl::save"
        "(FoamX::DictionaryWriter& dictWriter, bool dictEntry)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        word entryName = typeDescriptor_->name();
        string comment = typeDescriptor_->comment();
        string description = typeDescriptor_->description();

        if (dictEntry && debug::infoSwitch("FoamXwriteComments"))
        {
            if (comment.size() > 0)
            {
                dictWriter.writeComment(comment);
            }
            else
            {
                if (description.size() > 0)
                {
                    dictWriter.writeComment(description);
                }
            }
        }

        switch(typeDescriptor_->type())
        {
            case Type_FixedList:
            {
                log << "Saving FixedList" << endl;
            
                // Write entry name if required.
                if (dictEntry) dictWriter.writeKeyword(entryName);
            
                // Start list. Do not write the count.
                dictWriter.startFixedList();
            
                // Write values.
                label i = 0;
                for
                (
                    DLList<IDictionaryEntryImpl*>::iterator iter
                        = subElements_.begin();
                    iter != subElements_.end();
                    ++iter
                )
                {
                    if (i++ > 0)
                    {
                        dictWriter.writeString(" ");
                    }

                    iter()->save(dictWriter, false);
                }
            
                dictWriter.endFixedList();
            
                if (dictEntry)
                {
                    dictWriter.endEntry();
                }
            }
            break;

            case Type_List:
            {
                log << "Saving List" << endl;

                // Write entry name if required.
                if (dictEntry) dictWriter.writeKeyword(entryName);

                if (listTokenPtr_)
                {
                    dictWriter.writeValue(*listTokenPtr_);
                }
                else
                {
                    // Start list.
                    dictWriter.startList(subElements_.size());

                    // Write values.
                    for
                    (
                        DLList<IDictionaryEntryImpl*>::iterator iter = 
                            subElements_.begin();
                        iter != subElements_.end();
                        ++iter
                    )
                    {
                        if (iter()->typeDescriptor()->type() != Type_List)
                        {
                            dictWriter.writeString("\n");
                            dictWriter.indent();
                        }

                        iter()->save(dictWriter, false);
                    }

                    dictWriter.endList();
                }

                if (dictEntry) dictWriter.endEntry();
            }
            break;

            case Type_Dictionary:
            {
                log << "Saving Dictionary" << endl;

                // Start a sub-dictionary.
                if (dictEntry)
                {
                    dictWriter.startSubDict(entryName);
                }
                else
                {
                    dictWriter.startSubDict();
                }
            
                // Write sub-entries.
                for
                (
                    DLList<IDictionaryEntryImpl*>::iterator iter =
                        subElements_.begin();
                    iter != subElements_.end();
                    ++iter
                )
                {
                    iter()->save(dictWriter, true);
                }

                if (dictEntry) 
                {
                    dictWriter.endSubDict();
                }
                else
                {
                    dictWriter.endDict();
                }
            }
            break;

            case Type_Selection:
            {
                log << "Saving Selection" << endl;

                if (!subElements_.size())
                {
                    throw FoamXError
                    (
                        E_FAIL,
                        string("No types given for selection '")
                      + typeDescriptor_->path() + "' for dictionary\n"
                      + dictWriter.pathName(),
                        functionName,
                        __FILE__, __LINE__
                    );
                }

                // Write entry name if required.
                if (dictEntry) dictWriter.writeKeyword(entryName);

                label i = 0;

                // Write current selection.
                for
                (
                    DLList<IDictionaryEntryImpl*>::iterator iter = 
                        subElements_.begin();
                    iter != subElements_.end();
                    ++iter, ++i
                )
                {
                    if (i == selection_)
                    {
                        dictWriter.writeString
                        (
                            word(iter()->typeDescriptor_->name())
                        );

                        if
                        (
                            iter()->typeDescriptor_->isPrimitiveType()
                         || iter()->subElements_.size()
                         || iter()->listTokenPtr_
                        )
                        {
                            dictWriter.writeString(" ");
                            iter()->save(dictWriter, false);
                        }
                    }
                }
            
                if (dictEntry) dictWriter.endEntry();
            }
            break;

            case Type_Compound:
            {
                log << "Saving Compound" << endl;

                // Write entry name if required.
                if (dictEntry)
                {
                    if (subElements_.size())
                    {
                        dictWriter.writeKeyword(entryName);
                    }
                    else
                    {
                        dictWriter.writeKeywordOnly(entryName);
                    }
                }

                label i = 0;

                // Write values.
                for
                (
                    DLList<IDictionaryEntryImpl*>::iterator iter = 
                        subElements_.begin();
                    iter != subElements_.end();
                    ++iter
                )
                {
                    if (i++ > 0)
                    {
                        dictWriter.writeString(" ");
                    }

                    iter()->save(dictWriter, false);
                }
            
                if (dictEntry) dictWriter.endEntry();
            }
            break;

            case Type_Boolean:
            case Type_Label:
            case Type_Scalar:
            case Type_Char:
            case Type_Word:
            case Type_String:
            case Type_RootDir:
            case Type_RootAndCase:
            case Type_CaseName:
            case Type_HostName:
            case Type_File:
            case Type_Directory:
            case Type_Time:
            case Type_DimensionSet:
            {
                // If this entry is non-optional and a value hasn't been set,
                // throw our teddy right out of the pram.
                if (!typeDescriptor_->optional() && !value_.IsSet())
                {
                    throw FoamXError
                    (
                        E_FAIL,
                        string("Non-optional dictionary entry '")
                      + typeDescriptor_->path()
                      + "' : Value not specified for dictionary\n"
                      + dictWriter.pathName(),
                        functionName,
                        __FILE__, __LINE__
                    );
                }

                // If this entry is non-optional, or is optional and an entry
                // has been set, write it to the dictionary.
                if
                (
                    !typeDescriptor_->optional()
                 || (typeDescriptor_->optional() && value_.IsSet())
                )
                {
                    // Write entry name if required.
                    if (dictEntry) dictWriter.writeKeyword(entryName);

                    dictWriter.writeValue(value_);
                    if (dictEntry) dictWriter.endEntry();
                }
            }
            break;

            default:
                throw FoamXError
                (
                    E_INVALID_ARG,
                    "Invalid Type",
                    functionName,
                    __FILE__, __LINE__
                );
        }

        modified_ = false;
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// The DictionaryEntry object tree may have built with TypeDescriptor objects of
// type "field". These need to be replaced when we know the actual field type.
// The original type descriptor holds some type information regarding the
// original parameter (of type "field") such as the name, displayName and 
// Description.  Thus, we need to clone the field type descriptor and replace
// its name and description parameters with those originally specified.

void FoamX::IDictionaryEntryImpl::bindFieldType(ITypeDescriptor_ptr typeDesc)
{
    static const char* functionName =
        "FoamX::IDictionaryEntryImpl::bindFieldType"
        "(ITypeDescriptor_ptr typeDesc)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Check the type of the current type descriptor.
        if (typeDescriptor_->type() == Type_Field)
        {
            // Store original name, displayName and description parameters.
            CORBA::String_var name        = typeDescriptor_->name();
            CORBA::String_var displayName = typeDescriptor_->displayName();
            CORBA::String_var description = typeDescriptor_->description();

            // Clone the field type descriptor.
            ITypeDescriptorImpl* fieldTypeDescriptor = new ITypeDescriptorImpl
            (
                typeDesc
            );

            if (fieldTypeDescriptor == NULL)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Failed to create field TypeDescriptor object",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Copy over the name, displayName and description information.
            fieldTypeDescriptor->name(name);
            fieldTypeDescriptor->displayName(displayName);
            fieldTypeDescriptor->description(description);

            // Bind the given field type descriptor.
            bindType(fieldTypeDescriptor->_this());

            fieldTypeDescriptor->_remove_ref();
        }
        else if (typeDescriptor_->isCompoundType())
        {
            // Recurse over all sub-entries.
            for
            (
                DLList<IDictionaryEntryImpl*>::iterator iter =
                    subElements_.begin();
                iter != subElements_.end();
                ++iter
            )
            {
                iter()->bindFieldType(typeDesc);
            }
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IDictionaryEntryImpl::bindType
(
    FoamXServer::ITypeDescriptor_ptr typeDesc
)
{
    static const char* functionName =
        "FoamX::IDictionaryEntryImpl::bindType"
        "(FoamXServer::ITypeDescriptor_ptr typeDesc)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Reset type info.
        typeDescriptor_ = ITypeDescriptor::_nil();
        clearSubElements();

        // Duplicate and store the TypeDescriptor reference.
        typeDescriptor_ = ITypeDescriptor::_duplicate(typeDesc);

        log << "Constructing DictionaryEntry for "
            << CORBA::String_var(typeDescriptor_->path()) << endl;

        // See if this is a primitive type.
        if (typeDescriptor_->isPrimitiveType())
        {
            value_.setType(typeDescriptor_->type());
        }
        else
        {
            // Get list of the sub types from the type descriptor.
            // Auto release.
            TypeDescriptorList_var pSubTypes = typeDescriptor_->subTypes();

            switch(typeDescriptor_->type())
            {
                case Type_FixedList:
                {
                    log << "Default constructing FixedList"
                        << endl;
                
                    // Get size of vector.
                    label count = typeDescriptor_->numElements();
                
                    // A list type descriptor should only have one sub
                    // type - the type of the vector elements.
                    if (pSubTypes->length() != 1)
                    {
                        throw FoamXError
                        (
                            E_FAIL,
                            "Invalid number of SubTypes for vector, "
                            " expected 1, found "
                          + word
                            (
                                name(label(pSubTypes->length()))
                            ),
                            functionName,
                            __FILE__, __LINE__
                        );
                    }

                    // Add specified number of default elements.
                    for (int i = 0; i <count; i++)
                    {
                        IDictionaryEntryImpl* pSubEntry = new
                        IDictionaryEntryImpl(typeDescriptor_->elementType());

                        if (pSubEntry == NULL)
                        {
                            throw FoamXError
                            (
                                E_FAIL,
                                "Failed to create IDictionaryEntryImpl"
                                " object",
                                functionName,
                                __FILE__, __LINE__
                            );
                        }
                        subElements_.append(pSubEntry);
                    }
                }
                break;

                case Type_List:
                {
                    log << "Default constructing - List (with no elements)"
                        << endl;
                
                    // A list type descriptor should only have one sub
                    // type - the type of the list elements.
                    if (pSubTypes->length() != 1)
                    {
                        throw FoamXError
                        (
                            E_FAIL,
                            "Invalid number of SubTypes for list",
                            functionName,
                            __FILE__, __LINE__
                        );
                    }
                }
                break;

                case Type_Dictionary:
                {
                    log << "Default constructing Dictionary"
                        << endl;
                
                    // Loop over all sub types and create a default
                    // sub entry for each.
                    for (unsigned int i = 0; i <pSubTypes->length(); i++)
                    {
                        // Create default entry.
                        IDictionaryEntryImpl* pSubEntry = 
                            new IDictionaryEntryImpl(pSubTypes[i]);

                        if (pSubEntry == NULL)
                        {
                            throw FoamXError
                            (
                                E_FAIL,
                                "Failed to create IDictionaryEntryImpl"
                                " object",
                                functionName,
                                __FILE__, __LINE__
                            );
                        }
                        subElements_.append(pSubEntry);
                    }
                }
                break;

                case Type_Selection:
                {
                    log << "Default constructing Selection"
                        << endl;
                
                    // Loop over all sub types and create a default
                    // sub entry for each.
                    for (unsigned int i=0; i <pSubTypes->length(); i++)
                    {
                        // Create a default sub entry object and store its
                        // reference.
                        IDictionaryEntryImpl* pSubEntry = new
                        IDictionaryEntryImpl(pSubTypes[i]);

                        if (pSubEntry == NULL)
                        {
                            throw FoamXError
                            (
                                E_FAIL,
                                "Failed to create IDictionaryEntryImpl"
                                " object",
                                functionName,
                                __FILE__, __LINE__
                            );
                        }
                        subElements_.append(pSubEntry);
                    }
                }
                break;

                case Type_Compound:
                {
                    log << "Default constructing Compound"
                        << endl;
                
                    // Loop over all sub types and create a default
                    // sub entry for each.
                    for (unsigned int i=0; i <pSubTypes->length(); i++)
                    {
                        // Create a default sub entry object and store its
                        // reference.
                        IDictionaryEntryImpl* pSubEntry = new
                        IDictionaryEntryImpl(pSubTypes[i]);

                        if (pSubEntry == NULL)
                        {
                            throw FoamXError
                            (
                                E_FAIL,
                                "Failed to create IDictionaryEntryImpl"
                                " object",
                                functionName,
                                __FILE__, __LINE__
                            );
                        }
                        subElements_.append(pSubEntry);
                    }
                }
                break;

                case Type_Field:
                    log << "Default constructing Field" << endl;
                break;

                default:
                    throw FoamXError
                    (
                        E_INVALID_ARG,
                        "Invalid Type",
                        functionName,
                        __FILE__, __LINE__
                    );
            }
        }

        // Set this to defaultValue for any type
        if (typeDesc->hasDefaultValue())
        {
            IDictionaryEntry_var dval;
            typeDesc->getDefaultValue(dval.out());
            operator=(dval);
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IDictionaryEntryImpl::operator=
(
    const IDictionaryEntryImpl& dictEntry
)
{
    static const char* functionName =
        "FoamX::IDictionaryEntryImpl::operator=(const IDictionaryEntryImpl&)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        if (typeDescriptor_->type() != dictEntry.typeDescriptor_->type())
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Type of argument "
              + FoamXTypes::typeName(dictEntry.typeDescriptor_->type())
              + " does not match this type "
              + FoamXTypes::typeName(typeDescriptor_->type()),
                functionName,
                __FILE__, __LINE__
            );
        }

        value_ = dictEntry.value_;

        if (typeDescriptor_->type() == Type_List)
        {
            clearSubElements();

            for
            (
                DLList<IDictionaryEntryImpl*>::const_iterator
                    dictEntryIter = dictEntry.subElements_.begin();
                dictEntryIter != dictEntry.subElements_.end();
                ++dictEntryIter
            )
            {
                IDictionaryEntryImpl* pSubEntry = new IDictionaryEntryImpl
                (
                    *dictEntryIter()
                );

                if (pSubEntry == NULL)
                {
                    throw FoamXError
                    (
                        E_FAIL,
                        "Couldn't create IDictionaryEntryImpl object for "
                      + word(typeDescriptor_->name()),
                        functionName,
                        __FILE__, __LINE__
                    );
                }

                subElements_.append(pSubEntry);
            }

            if (listTokenPtr_)
            {
                delete listTokenPtr_;
                listTokenPtr_ = NULL;
            }

            if (dictEntry.listTokenPtr_)
            {
                listTokenPtr_ = new token(*dictEntry.listTokenPtr_);
            }
        }
        else
        {
            if (subElements_.size() != dictEntry.subElements_.size())
            {
                throw FoamXError
                (
                    E_INVALID_ARG,
                    "Number of sub-elements in argument "
                  + word(name(dictEntry.subElements_.size()))
                  + " does not match the number in this of "
                  + word(name(subElements_.size())),
                  functionName,
                  __FILE__, __LINE__
                );
            }

            DLList<IDictionaryEntryImpl*>::iterator iter =
                subElements_.begin();

            DLList<IDictionaryEntryImpl*>::const_iterator dictEntryIter =
                dictEntry.subElements_.begin();

            for
            (
                ;
                iter != subElements_.end()
             && dictEntryIter != dictEntry.subElements_.end();
                ++iter, ++dictEntryIter
            )
            {
                *iter() = *dictEntryIter();
            }

            if (typeDescriptor_->type() == Type_Selection)
            {
                selection_ = dictEntry.selection_;
            }
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IDictionaryEntryImpl::operator=(IDictionaryEntry_ptr dictEntryPtr)
{
    static const char* functionName =
        "FoamX::IDictionaryEntryImpl::operator=(IDictionaryEntry_ptr)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        if (typeDescriptor_->type() != dictEntryPtr->typeDescriptor()->type())
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Type of argument "
              + FoamXTypes::typeName(dictEntryPtr->typeDescriptor()->type())
              + " does not match this type "
              + FoamXTypes::typeName(typeDescriptor_->type()),
                functionName,
                __FILE__, __LINE__
            );
        }

        value_ = *dictEntryPtr->value();

        if (typeDescriptor_->type() == Type_List)
        {
            clearSubElements();

            FoamXServer::DictionaryEntryList* subElmtsPtr =
                dictEntryPtr->subElements();

            for (unsigned int i = 0; i < subElmtsPtr->length(); i++)
            {
                IDictionaryEntryImpl* pSubEntry = new IDictionaryEntryImpl
                (
                    (*subElmtsPtr)[i]->typeDescriptor()
                );

                if (pSubEntry == NULL)
                {
                    throw FoamXError
                    (
                        E_FAIL,
                        "Couldn't create IDictionaryEntryImpl object for "
                      + word((*subElmtsPtr)[i]->typeDescriptor()->name()),
                        functionName,
                        __FILE__, __LINE__
                    );
                }

                *pSubEntry = (*subElmtsPtr)[i];

                subElements_.append(pSubEntry);
            }
        }
        else
        {
            FoamXServer::DictionaryEntryList* subElmtsPtr =
                dictEntryPtr->subElements();

            if (subElements_.size() != label(subElmtsPtr->length()))
            {
                throw FoamXError
                (
                    E_INVALID_ARG,
                    "Number of sub-elements in argument "
                  + word(name(label(subElmtsPtr->length())))
                  + " does not match the number in this of "
                  + word(name(subElements_.size())),
                    functionName,
                    __FILE__, __LINE__
                );
            }

            label i = 0;
            for
            (
                DLList<IDictionaryEntryImpl*>::iterator iter =
                    subElements_.begin();
                iter != subElements_.end();
                ++iter
            )
            {
                *iter() = (*subElmtsPtr)[i++];
            }

            if (typeDescriptor_->type() == Type_Selection)
            {
                selection_ = dictEntryPtr->selection();
            }
        }
    }
    CATCH_ALL(functionName);
}


// ************************************************************************* //
