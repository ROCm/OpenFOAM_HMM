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
    Foam::IGeometricFieldDescriptorImpl

Description

SourceFiles
    IGeometricFieldImpl.C

\*---------------------------------------------------------------------------*/

#ifndef IGeometricFieldImpl_H
#define IGeometricFieldImpl_H

// Foam header files.
#include "fileName.H"
#include "dictionary.H"

// FoamX header files.
#include "FoamX.H"
#include "FoamXErrors.H"
#include "ObjRefHashTable.H"
#include "DictionaryWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace FoamX
{

// Forward declaration of classes.
class IDictionaryEntryImpl;

/*---------------------------------------------------------------------------*\
                   Class IGeometricFieldImpl Declaration
\*---------------------------------------------------------------------------*/

class IGeometricFieldImpl
:
    public POA_FoamXServer::CaseServer::IGeometricField,
    public PortableServer::RefCountServantBase
{
    // Private data

        //- Field name. No, really.
        Foam::word fieldName_;

        //- Reference to FoamProperties object. Used for patchField information.
        FoamXServer::CaseServer::IFoamProperties_var foamProperties_;

        //- Reference to FieldDescriptor object.
        FoamXServer::CaseServer::IGeometricFieldDescriptor_var fieldDescriptor_;

        //- Reference to internal field value object.
        IDictionaryEntryImpl* internalFieldValue_;

        //- Map of patch name to patch field types.
        Foam::HashTable<Foam::word> patchFieldNames_;

        //- Map of patch name to patch field parameter values.
        ObjRefHashTable<IDictionaryEntryImpl*> patchFields_;


public:

    // Constructors

        IGeometricFieldImpl
        (
            FoamXServer::CaseServer::IGeometricFieldDescriptor_ptr
                fieldDescriptor,
            FoamXServer::CaseServer::IFoamProperties_ptr foamProperties
        );


    // Destructor

            virtual ~IGeometricFieldImpl();


    // Member Functions
 
            virtual char* name();

            virtual void getInternalFieldValue
            (
                FoamXServer::IDictionaryEntry_out internalFieldValue
            );

            virtual void getPatchFieldParameters
            (
                const char* patchName,
                FoamXServer::IDictionaryEntry_out patchFieldParams
            );


        // Set patch field type for a particular patch.

            virtual void addPatch
            (
                const char* patchName,
                const char* patchFieldType
            );

            virtual void deletePatch(const char* patchName);

            virtual void setPatchFieldType
            (
                const char* patchName,
                const char* patchFieldType
            );


        // Recursive state functions.

            virtual CORBA::Boolean modified();


        // Persistence.

            void load(const Foam::dictionary& fieldDict);

            void save
            (
                DictionaryWriter& dictWriter,
                const Foam::wordList& patchNames
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FoamX

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
