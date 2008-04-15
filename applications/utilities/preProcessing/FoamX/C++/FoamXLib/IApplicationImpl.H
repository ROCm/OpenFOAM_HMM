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
    Foam::IApplicationImpl

Description

SourceFiles
    IApplicationImpl.C

\*---------------------------------------------------------------------------*/

#ifndef IApplicationImpl_H
#define IApplicationImpl_H

// Foam header files.
#include "fileName.H"
#include "dictionary.H"

// FoamX header files.
#include "FoamX.H"
#include "FoamXErrors.H"
#include "FoamXString.H"
#include "FoamXStringList.H"
#include "FoamXWordList.H"
#include "ObjRefDictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace FoamX
{

// Forward declaration of classes.
class IGeometricFieldDescriptorImpl;
class IPatchPhysicalTypeDescriptorImpl;
class ITypeDescriptorImpl;
class IPropertiesImpl;

/*---------------------------------------------------------------------------*\
                 Class IApplicationImpl Declaration
\*---------------------------------------------------------------------------*/

class IApplicationImpl
:
    public POA_FoamXServer::CaseServer::IApplication,
    public PortableServer::RefCountServantBase
{
    // Private data

        Foam::word name_;
        Foam::string description_;
        Foam::fileName category_;

        //- Path of the application class description
        Foam::fileName appCfgPath_;

        //- System or user class flag.
        bool systemClass_;

        //- Reference to the Foam properties class
        const IPropertiesImpl& foamSystemProperties_;

        ITypeDescriptorImpl* arguments_;
        FoamXWordList moduleList_;

        ObjRefDictionary<IGeometricFieldDescriptorImpl> fieldDescriptorMap_;

        ObjRefDictionary<IPatchPhysicalTypeDescriptorImpl>
            patchPhysicalTypeDescriptorMap_;

        ObjRefDictionary<ITypeDescriptorImpl> dictTypeDescriptorMap_;


public:

    // Constructors

        IApplicationImpl
        (
            const ApplicationDescriptor& appDesc,
            const IPropertiesImpl& foamSystemProperties
        );


    // Destructor

        virtual ~IApplicationImpl();


    // Member Functions

            virtual char* name();
            virtual void name(const char*);

            fileName appCfgPathName() const
            {
                return appCfgPath_/name_ + ".cfg";
            }

            virtual char* description();
            virtual void description(const char*);

            virtual char* category();
            virtual void category(const char*);

            virtual FoamXServer::StringList* modules();
            virtual void modules(const FoamXServer::StringList&);

            virtual CORBA::Boolean systemClass();

            virtual FoamXServer::StringList* fields();

            virtual void getField
            (
                const char* fieldKey,
                FoamXServer::CaseServer::IGeometricFieldDescriptor_out
                    fieldDescriptor
            );

            virtual void findField
            (
                const char* fieldName,
                FoamXServer::CaseServer::IGeometricFieldDescriptor_out
                    fieldDescriptor
            );

            virtual void addField
            (
                const char* fieldKey,
                FoamXServer::CaseServer::IGeometricFieldDescriptor_out
                    fieldDescriptor
            );

            virtual void deleteField(const char* fieldKey);

            virtual FoamXServer::StringList* patchPhysicalTypes();

            virtual void getPatchPhysicalType
            (
                const char* patchPhysicalTypeKey,
                FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_out
                    boundaryDescriptor
            );

            virtual void findPatchPhysicalType
            (
                const char* patchPhysicalTypeName,
                FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_out
                    boundaryDescriptor
            );

            virtual void addPatchPhysicalType
            (
                const char* patchPhysicalTypeKey,
                FoamXServer::CaseServer::IPatchPhysicalTypeDescriptor_out
                    boundaryDescriptor
            );

            virtual void deletePatchPhysicalType
            (
                const char* patchPhysicalTypeKey
            );

            virtual FoamXServer::StringList* dictionaries();

            virtual void getDictionary
            (
                const char* dictName,
                FoamXServer::ITypeDescriptor_out
                    dictTypeDescriptor
            );
    
            virtual void addDictionary
            (
                const char* dictName,
                FoamXServer::ITypeDescriptor_out
                    dictTypeDescriptor
            );

            virtual void deleteDictionary(const char* dictName);


        //- Validation.
        virtual void validate();

        //- Persistence.
        virtual void save();
 
        //- Load from application class definition dictionary.
        void load();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FoamX

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
