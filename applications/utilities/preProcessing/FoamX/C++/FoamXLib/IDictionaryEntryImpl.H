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
    Foam::IDictionaryEntryImpl

Description

SourceFiles
    IDictionaryEntryImpl.C

\*---------------------------------------------------------------------------*/

#ifndef IDictionaryEntryImpl_H
#define IDictionaryEntryImpl_H

// Foam header files.
#include "List.H"
#include "DLList.H"
#include "Istream.H"
#include "entry.H"

// FoamX header files.
#include "FoamX.H"
#include "FoamXErrors.H"
#include "FoamXAny.H"
#include "DictionaryWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace FoamX
{

// Forward declaration of classes.
class ITypeDescriptorImpl;

/*---------------------------------------------------------------------------*\
                  Class IDictionaryEntryImpl Declaration
\*---------------------------------------------------------------------------*/

class IDictionaryEntryImpl
:
    public POA_FoamXServer::IDictionaryEntry,
    public PortableServer::RefCountServantBase
{

protected:

    // Private data

        //- Type descriptor of this entry.
        FoamXServer::ITypeDescriptor_var typeDescriptor_;

        //- Modified flag.
        CORBA::Boolean modified_;

        //- Value of this entry for non-compound types.
        FoamXAny value_;

        //- Doubly linked list of sub-elements.
        Foam::DLList<IDictionaryEntryImpl*> subElements_;

        //- Pointer to hold the list token for large lists
        Foam::token* listTokenPtr_;

        //- Current selection
        label selection_;


    // Private member functions.

        void clearSubElements();

        template<class Type, class CorbaType>
        void expandPrimitiveList();

        template<class Form>
        void expandFixedListList();

        void expandList();

        void bindType
        (
            FoamXServer::ITypeDescriptor_ptr typeDesc
        );


public:

    // Constructors

        //- Construct default for given type descriptor
        IDictionaryEntryImpl(FoamXServer::ITypeDescriptor_ptr typeDesc);

        //- Construct as copy
        IDictionaryEntryImpl(IDictionaryEntryImpl&);


    // Destructor

        virtual ~IDictionaryEntryImpl();


    // Member Functions
 
        virtual FoamXServer::ITypeDescriptor_ptr typeDescriptor();


        // Non-compound type functions.

            virtual FoamXServer::FoamXAny* value();
            virtual void value(const FoamXServer::FoamXAny& newValue);
            virtual void setValue(const FoamXServer::FoamXAny& newValue);


        // Compound type functions.

            const Foam::DLList<IDictionaryEntryImpl*>& subElementPtrs() const
            {
                return subElements_;
            }

            virtual FoamXServer::DictionaryEntryList* subElements();
            virtual CORBA::Long nSubElements();
            virtual CORBA::Boolean packedList();

        // Dictionary type functions.

//            virtual CORBA::Boolean found(const char* entryName) const;
//            virtual FoamXServer::DictionaryEntry* lookup
//            (
//                const char* entryName
//            ) const;
//            virtual FoamXServer::DictionaryEntry* subDict
//            (
//                const char* entryName
//            ) const;


        // Selection type functions.

            virtual CORBA::Long selection();
            virtual void selection(CORBA::Long);


        // List type functions.

            virtual void addElement
            (
                FoamXServer::IDictionaryEntry_out subEntry
            );

            virtual void removeElement
            (
                FoamXServer::IDictionaryEntry_ptr subEntry
            );


        // Recursive state functions.

            virtual void validate();
            virtual CORBA::Boolean modified();


        // Persistence.

            void load(const Foam::entry&);
            void load(const Foam::dictionary&, bool allowNonOptional = false);
            virtual void load(Foam::Istream&);

            //- Save top-level dictionary.
            virtual void save();

            void save
            (
                DictionaryWriter& dictWriter,
                bool dictEntry = false
            );


        //- Method to allow field type descriptors to be late-bound.
        void bindFieldType
        (
            FoamXServer::ITypeDescriptor_ptr typeDesc
        );


    // Member operators

        void operator=(const IDictionaryEntryImpl&);
        void operator=(IDictionaryEntry_ptr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FoamX

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
