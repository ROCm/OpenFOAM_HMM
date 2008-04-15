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
    Foam::ICaseServerImpl

Description

SourceFiles
    ICaseServerImpl.C

\*---------------------------------------------------------------------------*/

#ifndef ICaseServerImpl_H
#define ICaseServerImpl_H

// Foam header files.
#include "word.H"
#include "stringList.H"
#include "fileName.H"
#include "Dictionary.H"
#include "Time.H"

// FoamX header files.
#include "FoamX.H"
#include "ProcessControl.H"
#include "ObjRefHashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace FoamX
{

// Forward declaration of classes.
class IGeometricFieldImpl;
class RootDictionary;
class PatchProperties;
class Orb;

/*---------------------------------------------------------------------------*\
                         Class ICaseServerImpl Declaration
\*---------------------------------------------------------------------------*/

class ICaseServerImpl
:
    public POA_FoamXServer::CaseServer::ICaseServer,
    public PortableServer::RefCountServantBase
{
    // Private data

        //- Case root directory.
        Foam::fileName rootDir_;

        //- Case name.
        Foam::fileName caseName_;

        //- Context name for name server registration.
        Foam::fileName hostContext_; 

        //- Context name for name server registration.
        Foam::fileName userContext_; 

        //- Object name for name server registration.
        Foam::fileName objectName_; 

        //- Case dictionary filename.
        Foam::fileName caseDictName_;

        //- Name of the application class for this case.
        Foam::word appName_;

        //- Reference to the ORB
        Orb& orb_;

        //- Current database (only in 'open' mode)
        Foam::Time* dbPtr_;

        //- Reference to the case browser
        FoamXServer::CaseBrowser::ICaseBrowser_var caseBrowser_;

        //- Reference to the FoamProperties object.
        FoamXServer::CaseServer::IFoamProperties_var foamProperties_;

        //- Reference to the Application object for this case.
        FoamXServer::CaseServer::IApplication_var app_;

        //- Managed flag.
        bool managed_;

        //- Process controls
        FoamX::ProcessControl procControl_;

        //- Pid of running case.
        pid_t pid_;

        //- Managed Dictionary of patch name to PatchProperties objects.
        Foam::Dictionary<PatchProperties> patchMap_;

        //- Map of field name to FieldValues objects.
        ObjRefHashTable<IGeometricFieldImpl*> fieldValueMap_;

        //- Map of dictionary name to Dictionary objects.
        ObjRefHashTable<RootDictionary*> dictionaryMap_;


    // Private member functions.

        //- Access database
        Foam::Time& time();

        //- Save (default) controlDict
        void saveControlDict();

        void writePatchData();
        void addPatch(PatchProperties* patchProps);
        bool fieldsMatchPatchPhysicalType(const Foam::word& patchPhysicalType);


public:

    // Constructors

        ICaseServerImpl
        (
            Orb& orb,
            const Foam::fileName& rootDir,
            const Foam::fileName& caseName,
            const Foam::word& mode,
            const Foam::word& app
        );


    // Destructor

        virtual ~ICaseServerImpl();


    // Member Functions
 
        // Implemented ICaseServerImpl functions.

            virtual CORBA::Boolean managed();
            virtual void managed(CORBA::Boolean);
            virtual char* caseRoot();
            virtual char* caseName();


        //- Application class for this case (read only).
        virtual FoamXServer::CaseServer::IApplication_ptr application();

        //- Foam system properties (read only).
        virtual FoamXServer::CaseServer::IFoamProperties_ptr foamProperties();

        //- Time-steps.
        virtual FoamXServer::StringList* availableTimeSteps();
        virtual char* getTime();
        virtual void setTime(const char*, const CORBA::Long);

        // Mesh.

            virtual CORBA::Boolean meshDefined();
            virtual void readMesh();
            virtual void importMesh
            (
                const char* hostname,
                const char* rootDir,
                const char* caseName
            );


        //- Field values.
        virtual void getFieldValues
        (
            const char* fieldName,
            FoamXServer::CaseServer::IGeometricField_out fieldValues
        );


        // Patches.

            virtual FoamXServer::StringList* patchNames();

            virtual void addPatch
            (
                const char* patchName,
                const char* patchPhysicalType
            );

            virtual void deletePatch(const char* patchName);

            virtual void deleteAllPatches();


        // Patch boundary conditions.

            virtual void setPatchPhysicalType
            (
                const char* patchName,
                const char* patchPhysicalType
            );

            virtual void getPatchPhysicalType
            (
                const char* patchName,
                CORBA::String_out patchPhysicalType
            );


        // Dictionaries.

            virtual void getDictionary
            (
                const char* dictionaryName,
                CORBA::Boolean,
                FoamXServer::IDictionaryEntry_out dictRoot
            );

            virtual CORBA::Long fileModificationDate(const char* dictName);

            virtual void readFile
            (
                const char* name,
                CORBA::String_out contents
            );

            virtual void writeFile
            (
                const char* name,
                const char* contents
            );

        //- Calculation control.
        virtual CORBA::Long runCase(const char* arguments);

        virtual void killCase();

        //- Validation.
        virtual void validate();

        //- Has it been modified?
        virtual CORBA::Boolean modified();

        //- Persistence.
        virtual void save();

        //- Lifetime management.
        virtual void close();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FoamX

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
