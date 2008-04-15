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
    Foam::IPatchPhysicalTypeDescriptorImpl

Description

SourceFiles
    IPatchPhysicalTypeDescriptorImpl.C

\*---------------------------------------------------------------------------*/

#ifndef IPatchPhysicalTypeDescriptorImpl_H
#define IPatchPhysicalTypeDescriptorImpl_H

// Foam header files.
#include "IDLList.H"
#include "dictionary.H"

// FoamX header files.
#include "FoamX.H"
#include "DictionaryWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace FoamX
{

/*---------------------------------------------------------------------------*\
              Class IPatchPhysicalTypeDescriptorImpl Declaration
\*---------------------------------------------------------------------------*/

class IPatchPhysicalTypeDescriptorImpl
:
    public POA_FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor,
    public PortableServer::RefCountServantBase,
    public IDLList<IPatchPhysicalTypeDescriptorImpl>::link
{
    // Private data

        Foam::word name_;
        Foam::string displayName_;
        Foam::string description_;
        Foam::word patchType_;
        Foam::word parentType_;
        FoamXServer::StringPairList patchFieldTypes_;


    // Private member functions

        //- Search back through the parentTypes to find the patchPhysicalType
        //  for which the patchType is defined
        bool findParentPatchType
        (
            word& patchType,
            const word& parentType,
            const dictionary& patchPhysicalTypesDict
        ) const;

        //- Search back through the parentTypes to find the patchPhysicalType
        //  for which this boundary condition for this field is defined
        bool findPatchFieldType
        (
            word& patchFieldType,
            const word& parentType,
            const dictionary& patchPhysicalTypesDict,
            const dictionary& patchFieldPhysicalTypes
        ) const;


public:

    // Constructors

        IPatchPhysicalTypeDescriptorImpl
        (
            const word& patchPhysicalTypeName
        );

        IPatchPhysicalTypeDescriptorImpl
        (
            const word& patchPhysicalTypeName,
            const Foam::wordList& fieldList
        );

        IPatchPhysicalTypeDescriptorImpl
        (
            const word& patchPhysicalTypeName,
            const Foam::dictionary& patchPhysicalTypes,
            const dictionary& patchPhysicalTypesDict,
            const Foam::dictionary& patchFieldsPhysicalTypesDict,
            const Foam::wordList& fieldList
        );


    // Destructor

        virtual ~IPatchPhysicalTypeDescriptorImpl();


    // Member Functions

        // Foam access functions

            //- Keyword used by Dictionary
            const word& keyword() const
            {
                return name_;
            }


        // FoamServer access functions

            virtual char* name();
            virtual void  name(const char*);
            virtual char* displayName();
            virtual void  displayName(const char*);
            virtual char* description();
            virtual void description(const char*);
            virtual char* patchType();
            virtual void patchType(const char*);
            virtual char* parentType();
            virtual void parentType(const char*);
            virtual FoamXServer::StringPairList* patchFieldTypes();
            virtual void patchFieldTypes(const FoamXServer::StringPairList&);

            //- Persist to file.
            void save(DictionaryWriter& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FoamX

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
