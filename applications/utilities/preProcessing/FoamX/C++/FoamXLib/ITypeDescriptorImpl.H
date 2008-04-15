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

Class
    Foam::ITypeDescriptorImpl

Description

SourceFiles
    ITypeDescriptorImpl.C

\*---------------------------------------------------------------------------*/

#ifndef ITypeDescriptorImpl_H
#define ITypeDescriptorImpl_H

// Foam header files.
#include "IDLList.H"
#include "DLList.H"
#include "dictionary.H"

// FoamX header files.
#include "FoamX.H"
#include "FoamXErrors.H"
#include "FoamXAny.H"
#include "IDictionaryEntryImpl.H"
#include "FoamXWordList.H"
#include "FoamXTypes.H"
#include "DictionaryWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace FoamX
{

/*---------------------------------------------------------------------------*\
                  Class ITypeDescriptorImpl Declaration
\*---------------------------------------------------------------------------*/

class ITypeDescriptorImpl
:
    public POA_FoamXServer::ITypeDescriptor,
    public PortableServer::RefCountServantBase,
    public IDLList<ITypeDescriptorImpl>::link
{
    // Private data

        //- The type of this parameter.
        FoamXServer::FoamXType type_;

        //- Foam parameter name.
        Foam::word name_;

        //- Path of this parameter relative to root dictionary.
        Foam::string path_;

        //- Short display name.
        Foam::string displayName_;

        //- Long descriptive name.
        Foam::string description_;

        //- Dictionary comment.
        Foam::string comment_;

        //- GUI category for this parameter.
        Foam::string category_;

        //- URL to relevant help file.
        Foam::string helpURL_;

        //- URL to a nice icon for this parameter.
        Foam::string iconURL_;

        //- Flag indicating whether this parameter is optional.
        //  Defaults to false.
        bool optional_;

        //- Flag indicating whether this parameter is visible within FoamX.
        //  Defaults to true.
        bool visible_;

        //- Flag indicating whether this parameter is editable within FoamX.
        //  Defaults to true.
        bool editable_;

        //- Minimum value for this parameter.
        FoamXAny minValue_; 

        //- Maximum value for this parameter.
        FoamXAny maxValue_;

        //- Name of dictionary to use as a lookup table.
        FoamXAnyList valueList_;

        //- List of permissible values for this parameter.
        Foam::string lookupDict_; 

        //- Top-level dictionary location relative to case root.
        Foam::string dictionaryPath_;

        //- Number of vector space elements.
        Foam::label numElements_;

        //- Labels of vector space elements (eg, x, y, z).
        FoamXWordList elementLabels_;

        //- Doubly linked list of sub-type objects. Used to maintain the order
        //  of sub-types.
        Foam::DLList<ITypeDescriptorImpl*> subTypes_;

        //- Default value for this parameter.
        IDictionaryEntryImpl* default_;


    // Private member functions.
    
        bool setType(const Foam::word& type);

        void resolveType
        (
            Foam::word& type,
            const Foam::dictionary& typeDict,
            const Foam::dictionary& foamTypesDict
        );

        void resolveTypeDict
        (
            Foam::word& type,
            const Foam::dictionary& typeDict,
            const Foam::dictionary& foamTypesDict
        );

        void readType
        (
            Foam::word& type,
            const Foam::dictionary& typeDict,
            const Foam::dictionary& foamTypesDict
        );

        void addEntries
        (
            const Foam::dictionary& typeDict,
            const Foam::dictionary& foamTypesDict
        );

        void addDictionaryEntries
        (
            const Foam::dictionary& typeDict,
            const Foam::dictionary& foamTypesDict
        );

        void addCompoundEntries
        (
            const Foam::dictionary& typeDict,
            const Foam::dictionary& foamTypesDict
        );

        void addElementType
        (
            const Foam::dictionary& typeDict,
            const Foam::dictionary& foamTypesDict
        );

        void readOptionalData
        (
            const Foam::dictionary& typeDict
        );


public:

    // Constructors

        ITypeDescriptorImpl
        (
            const Foam::word& name,
            const FoamXServer::FoamXType& type,
            const Foam::string& parentPath
        );

        ITypeDescriptorImpl
        (
            const Foam::word& name,
            const Foam::string& parentPath,
            const Foam::dictionary& typeDict,
            const Foam::dictionary& foamTypesDict
        );

        ITypeDescriptorImpl
        (
            const Foam::word& name,
            const Foam::string& parentPath,
            const Foam::entry& typeEntry,
            const Foam::dictionary& foamTypesDict
        );

        ITypeDescriptorImpl
        (
            FoamXServer::ITypeDescriptor_ptr typeDesc
        );


    // Destructor

        virtual ~ITypeDescriptorImpl();


    // Member Functions
 
        // Foam access functions

            //- Keyword used by Dictionary
            const word& keyword() const
            {
                return name_;
            }


        // FoamServer access functions

            virtual FoamXServer::FoamXType type();
            virtual void type(FoamXServer::FoamXType);
            virtual char* path();
            virtual void path(const char*);
            virtual char* name();
            virtual void name(const char*);
            virtual char* displayName();
            virtual void displayName(const char*);
            virtual char* description();
            virtual void description(const char*);
            virtual char* comment();
            virtual void comment(const char*);
            virtual char* category();
            virtual void category(const char*);
            virtual char* helpURL();
            virtual void helpURL(const char*);
            virtual char* iconURL();
            virtual void iconURL(const char*);
            virtual CORBA::Boolean isPrimitiveType();
            virtual CORBA::Boolean isCompoundType();
            virtual CORBA::Boolean optional();
            virtual void optional(CORBA::Boolean);
            virtual CORBA::Boolean visible();
            virtual void visible(CORBA::Boolean);
            virtual CORBA::Boolean editable();
            virtual void editable(CORBA::Boolean);


        // Properties for non-compound types.

            virtual FoamXServer::FoamXAny* minValue();
            virtual void minValue(const FoamXServer::FoamXAny&);
            virtual FoamXServer::FoamXAny* maxValue();
            virtual void maxValue(const FoamXServer::FoamXAny&);
            virtual FoamXServer::FoamXAnyList* valueList();
            virtual void valueList(const FoamXServer::FoamXAnyList&);
            virtual char* lookupDict();
            virtual void lookupDict(const char*);


        // Properties for compound types.

            virtual char* dictionaryPath();
            virtual void dictionaryPath(const char*);
            virtual CORBA::Long numElements();
            virtual void numElements(CORBA::Long);
            virtual FoamXServer::StringList* elementLabels();
            virtual void elementLabels(const FoamXServer::StringList&);
            virtual FoamXServer::TypeDescriptorList* subTypes();
            virtual FoamXServer::ITypeDescriptor_ptr elementType();


        // Default handling

            virtual CORBA::Boolean hasDefaultValue();

            virtual void getDefaultValue
            (
                FoamXServer::IDictionaryEntry_out defaultValue
            );


        // Compound sub-type methods.

            virtual void addSubType
            (
                FoamXServer::FoamXType type,
                FoamXServer::ITypeDescriptor_out subEntry
            );

            virtual void removeSubType
            (
                FoamXServer::ITypeDescriptor_ptr subEntry
            );


        //- Validation.
        virtual void validate();


    // Persistence.

        void load
        (
            const Foam::dictionary& typeDict,
            const Foam::dictionary& foamTypesDict
        );

        void save
        (
            DictionaryWriter& dictWriter,
            bool topLevelDict
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FoamX

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
