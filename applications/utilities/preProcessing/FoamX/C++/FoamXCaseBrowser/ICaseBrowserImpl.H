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
    ICaseBrowserImpl.C

\*---------------------------------------------------------------------------*/

#ifndef ICaseBrowserImpl_H
#define ICaseBrowserImpl_H

// Foam header files.
#include "fileName.H"
#include "HashPtrTable.H"
#include "long.H"
#include "stringList.H"

// FoamX header files.
#include "FoamX.H"
#include "ProcessControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    inline bool operator==
    (
        const FoamXServer::JobID& jid1,
        const FoamXServer::JobID& jid2
    )
    {
        return
        (
            jid1.processID == jid2.processID
         && Foam::word(jid1.hostName) == Foam::word(jid2.hostName)
        );
    }

    inline Ostream& operator<<
    (
        Ostream& os,
        const FoamXServer::JobID& jid
    )
    {
        os << Foam::word(jid.hostName) << '.' << jid.processID;
        
        return os;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace FoamX
{

// Forward declaration of classes.
class IPropertiesImpl;
class Orb;

/*---------------------------------------------------------------------------*\
                           Class ICaseBrowserImpl Declaration
\*---------------------------------------------------------------------------*/

class ICaseBrowserImpl
:
    public POA_FoamXServer::CaseBrowser::ICaseBrowser,
    public PortableServer::RefCountServantBase
{
    // Private data

        //- Context name for name server registration.
        Foam::fileName hostContext_;

        //- Object name for name server registration.
        Foam::fileName objectName_;

        //- Reference to the ORB
        Orb& orb_;

        //- Argument list
        Foam::stringList args_;

        //- Reference to the host browser
        FoamXServer::HostBrowser::IHostBrowser_var hostBrowser_;

        //- Foam system properties for this server.
        IPropertiesImpl* foamProperties_;

        //- Case list.
        Foam::HashPtrTable<FoamXServer::CaseDescriptor, Foam::string> cases_;

        class hashJobID
        {

        public:

            hashJobID()
            {}

            Foam::label operator()
            (
                const FoamXServer::JobID& key,
                const Foam::label tableSize
            ) const
            {
                return Foam::mag(Foam::label(key.processID%tableSize));
            }
        };

        typedef Foam::HashPtrTable
            <FoamXServer::JobDescriptor, FoamXServer::JobID, hashJobID>
            JobHashTable;

        //- Running jobs list.
        JobHashTable runningJobs_;

        //- Finished jobs list.
        JobHashTable finishedJobs_;

        //- Process controls
        FoamX::ProcessControl procControl_;


    // Private member functions.

        // Return the case server key for the given root and case
        static Foam::fileName caseServerKey
        (
            const char* rootDir,
            const char* caseName
        );

        void getCaseServerReference
        (
            const Foam::fileName& lockFile,
            FoamXServer::CaseServer::ICaseServer_out caseObj
        );

        static FoamXServer::CaseDescriptor* caseDescriptor
        (
            const Foam::fileName& rootDir,
            const Foam::fileName& rawRootDir,
            const Foam::fileName& caseFileName
        );

        //- Get dictionary entry or default value.
        static Foam::label getEntry
        (
            const Foam::dictionary&,
            const Foam::word& key,
            const Foam::label defaultValue
        );

        //- Get dictionary entry or default value.
        static Foam::string getEntry
        (
            const Foam::dictionary&,
            const Foam::word& key,
            const Foam::string& defaultValue
        );

        static void readJobs
        (
            JobHashTable& jobsList,
            const Foam::fileName& dir,
            const bool throwOnError
        );

public:

    // Constructors

        ICaseBrowserImpl(Orb& orb, const Foam::stringList& args);


    // Destructor

        virtual ~ICaseBrowserImpl();


    // Member Functions

        // Implemented ICaseBrowser methods.

            //- Foam system properties (user data editable).
            virtual FoamXServer::CaseServer::IFoamProperties_ptr
                foamProperties();

            //- Managed cases list.
            virtual FoamXServer::CaseDescriptorList* cases();

            //- Running jobs list.
            virtual FoamXServer::JobDescriptorList* runningJobs();

            //- Finished jobs list.
            virtual FoamXServer::JobDescriptorList* finishedJobs();

            virtual void getEnv
            (
                const char* envName,
                CORBA::String_out enval
            );

            virtual void getHostName(CORBA::String_out hostName);

            virtual void getUserName(CORBA::String_out userName);

            virtual CORBA::Long fileModificationDate(const char* fName);

            virtual void readFile
            (
                const char* fName,
                CORBA::String_out contents
            );

            virtual void writeFile
            (
                const char* fName,
                const char* contents
            );

            virtual CORBA::Long invokeUtility
            (
                const char* hostName,
                const char* utilityName,
                const FoamXServer::StringList& arguments,
                const char* logName,
                CORBA::Boolean background
            );

            virtual void openCase
            (
                const FoamXServer::CaseDescriptor& caseDesc
            );

            virtual void openCasePost
            (
                const FoamXServer::CaseDescriptor& caseDesc,
                const CORBA::Long nProcs
            );

            virtual void newCase
            (
                const char* rootDir,
                const char* caseName,
                const char* app
            );

            virtual void importCase
            (
                const char* rootDir,
                const char* caseName,
                const char* app
            );

            virtual void deleteCase
            (
                const FoamXServer::CaseDescriptor& caseDesc
            );

            virtual void cloneCase
            (
                const FoamXServer::CaseDescriptor& caseDesc,
                const char* newRootDir,
                const char* newCaseName,
                const char* newAppClassName,
                const char* timeSel
            );

            virtual CORBA::Boolean getCaseServerReference
            (
                const char* rootDir,
                const char* caseName,
                FoamXServer::CaseServer::ICaseServer_out caseObj
            );

            virtual CORBA::Boolean getCasePostServerReference
            (
                const char* rootDir,
                const char* caseName,
                CORBA::Long nProcs,
                FoamXServer::CasePostServer::ICasePostServer_out caseObj
            );

            virtual CORBA::Boolean caseLocked
            (
                const FoamXServer::CaseDescriptor& caseDesc
            );

            virtual void unlockCase
            (
                const char* rootDir,
                const char* caseName
            );

            virtual void unlockCaseDescriptor
            (
                const FoamXServer::CaseDescriptor& caseDesc
            );

            virtual void addCase
            (
                const char* rootDir,
                const char* rawRootDir,
                const char* caseName,
                const char* appName
            );

            virtual void caseOpen
            (
                const char* rootDir,
                const char* caseName
            );

            virtual CORBA::Boolean isCaseInError
            (
                const FoamXServer::CaseDescriptor& caseDesc
            );

            virtual void caseIsInError
            (
                const FoamXServer::CaseDescriptor& caseDesc
            );

            virtual void refreshCaseList();

            virtual void addToCaseList(const char* rootDir);


        // Process control

            virtual void refreshJobsLists();
            virtual void checkRunningJobs();
            virtual void purgeRunningJobs();
            virtual void purgeFinishedJob(const FoamXServer::JobID& jobID);
            virtual void purgeFinishedJobs(CORBA::Long nDays);
            virtual void kill(const FoamXServer::JobID& jobID);
            virtual void suspend(const FoamXServer::JobID& jobID);
            virtual void cont(const FoamXServer::JobID& jobID);
            virtual void end
            (
                const FoamXServer::JobID& jobID,
                const char* rootDir,
                const char* caseName,
                const CORBA::Boolean now                
            );
            virtual void setStatus
            (
                const FoamXServer::JobID& jobID,
                FoamXServer::JobStatus jobStatus
            );


       // IO

            virtual void validate();
            virtual void save();
            virtual void close();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FoamX

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
