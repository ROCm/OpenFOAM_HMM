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
#include "long.H"
#include "dictionary.H"
#include "IFstream.H"
#include "OFstream.H"
#include "wordList.H"
#include "stringList.H"
#include "fileNameList.H"
#include "clock.H"
#include "OSspecific.H"
#include "Pstream.H"
#include "HashSet.H"

// FoamX header files.
#include "ICaseBrowserImpl.H"
#include "IPropertiesImpl.H"
#include "LogEntry.H"
#include "NameServer.H"
#include "DictionaryWriter.H"
#include "Orb.H"
#include "Paths.H"
#include "instantList.H"
#include "Time.H"
#include "Time.H"

// Namespaces
#include "FoamXNameSpaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::ICaseBrowserImpl::ICaseBrowserImpl
(
    Orb& orb,
    const stringList& args
)
:
    hostContext_(hostName()),
    objectName_(hostName()/"FoamXCaseBrowser"),
    orb_(orb),
    args_(args),
    hostBrowser_(NULL),
    foamProperties_(NULL),
    procControl_(dotFoam("apps/FoamX/FoamX.cfg"))
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::ICaseBrowserImpl(Orb& orb)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Create the Foam system properties object.
        // Allow application class authoring.
        foamProperties_ = new IPropertiesImpl(false);
        if (foamProperties_ == NULL)
        {
            throw FoamXError
            (
                E_FAIL,
                "Failed to create Foam system properties object",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Connect to the naming service.
        NameServer fxNameServer(orb_.orbPtr_);

        // First get reference to hostBrowser
        // Host browser should have registered itself with the name server
        // under a "FoamXHostBrowser" key.
        fileName hostBrowserKey = "FoamXHostBrowser";

        if (fxNameServer.isObjectBound(hostBrowserKey))
        {
            // Return server reference.
            hostBrowser_ = 
                fxNameServer.resolve<FoamXServer::HostBrowser::IHostBrowser>
                (hostBrowserKey);
        }
        else
        {
            throw FoamXError
            (
                E_FAIL,
                "Hostbrowser '" + hostBrowserKey
              + "' not found on name server. "
                "getServerReference call timed out",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Register myself with the naming service.

        // Do not fail if context already exists.
        fxNameServer.createContext(hostContext_, false);

        // Bind object under key hostName/CaseBrowser.
        FoamXServer::CaseBrowser::ICaseBrowser_var browserRef = _this();

        // Throw exception if already bound.
        fxNameServer.bindObject(objectName_, browserRef, false);

        // Decrement the reference count of the object implementation, so
        // that it will be properly cleaned up when the POA has determined
        // that it is no longer needed.
        _remove_ref();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::ICaseBrowserImpl::~ICaseBrowserImpl()
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::~ICaseBrowserImpl()";

    LogEntry log(functionName, __FILE__, __LINE__);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::CaseServer::IFoamProperties_ptr
FoamX::ICaseBrowserImpl::foamProperties()
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::foamProperties()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Return a reference to the Foam System Properties object.
    IFoamProperties_ptr pFoamProperties = IFoamProperties::_nil();

    if (foamProperties_ != NULL)
    {
        pFoamProperties = foamProperties_->_this();
    }

    return pFoamProperties;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::CaseDescriptorList* FoamX::ICaseBrowserImpl::cases()
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::cases()";

    LogEntry log(functionName, __FILE__, __LINE__);

    CaseDescriptorList* caseListPtr = new CaseDescriptorList();
    caseListPtr->length(cases_.size());

    label i = 0;
    for
    (
        Foam::HashPtrTable<FoamXServer::CaseDescriptor, string>::iterator
            iter = cases_.begin();
        iter != cases_.end();
        ++iter
    )
    {
        (*caseListPtr)[i++] = *iter();
    }

    return caseListPtr;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::JobDescriptorList* FoamX::ICaseBrowserImpl::runningJobs()
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::runningJobs()";

    LogEntry log(functionName, __FILE__, __LINE__);

    JobDescriptorList* jobListPtr = new JobDescriptorList();
    jobListPtr->length(runningJobs_.size());

    label i = 0;
    for
    (
        JobHashTable::iterator iter = runningJobs_.begin();
        iter != runningJobs_.end();
        ++iter
    )
    {
        (*jobListPtr)[i++] = *iter();
    }

    return jobListPtr;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::JobDescriptorList* FoamX::ICaseBrowserImpl::finishedJobs()
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::finishedJobs()";

    LogEntry log(functionName, __FILE__, __LINE__);

    JobDescriptorList* jobListPtr = new JobDescriptorList();
    jobListPtr->length(finishedJobs_.size());

    label i = 0;
    for
    (
        JobHashTable::iterator
            iter = finishedJobs_.begin();
        iter != finishedJobs_.end();
        ++iter
    )
    {
        (*jobListPtr)[i++] = *iter();
    }

    return jobListPtr;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::getEnv
(
    const char* envName,
    CORBA::String_out enval
)
{
    enval = CORBA::string_dup(Foam::getEnv(envName).c_str());
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::getHostName
(
    CORBA::String_out hostName
)
{
    hostName = CORBA::string_dup(Foam::hostName().c_str());
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::getUserName
(
    CORBA::String_out userName
)
{
    userName = CORBA::string_dup(Foam::userName().c_str());
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

CORBA::Long FoamX::ICaseBrowserImpl::fileModificationDate
(
    const char* fName
)
{
    if (exists(fName))
    {
        return label(lastModified(fName));
    }
    else
    {
        //mj. Ugly. Let's hope not many dates with time_t==0
        return 0;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::readFile
(
    const char* fName,
    CORBA::String_out contents
)
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::readFile"
        "(const char* fName, CORBA::String_out contents)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        fileName ffName = fName;

        log << "File = " << ffName << endl;

        // If file exists, return its contents.
        if (exists(ffName))
        {
            IFstream ifFile(ffName);

            if (ifFile)
            {
                OStringStream fileString;

                char ch;
                while (ifFile.get(ch))
                {
                    fileString << ch;
                }

                // Must duplicate string to return.
                contents = CORBA::string_dup(fileString.str().c_str());
            }
            else
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Failed to read file '" + ffName + "'",
                    functionName,
                    __FILE__, __LINE__
                );
            }
        }
        else
        {
            throw FoamXError
            (
                E_FAIL,
                "Failed to open file '" + ffName + "'",
                functionName,
                __FILE__, __LINE__
            );
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::writeFile
(
    const char* fName,
    const char* contents
)
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::writeFile"
        "(const char* fName, const char* contents)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        fileName ffName = fName;

        log << "File = " << ffName << endl;

        // Check if file already exists, and if so move to .bak
        if (file(ffName))
        {
            mv(ffName, ffName + ".bak");
        }

        // Make sure path exists.
        if (!dir(ffName.path()) && !mkDir(ffName.path()))
        {
            throw FoamXError
            (
                E_FAIL,
                "Failed to create directory '" + ffName.path() + "'",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Write file contents.
        OFstream ofFile(ffName);
        if (ofFile)
        {
            ofFile << contents;
        }
        else
        {
            throw FoamXError
            (
                E_FAIL,
                "Failed to write file '" + ffName + "'",
                functionName,
                __FILE__, __LINE__
            );
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

CORBA::Long FoamX::ICaseBrowserImpl::invokeUtility
(
    const char* hostName,
    const char* utilityName,
    const FoamXServer::StringList& arguments,
    const char* logName,
    CORBA::Boolean background
)
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::invokeUtility"
        "(const char* hostName, const char* utilityName, "
        "const FoamXServer::StringList& arguments, const char* logName"
        "CORBA::Boolean background)";

    LogEntry log(functionName, __FILE__, __LINE__);

    pid_t pid = -1;

    try
    {
        log << "Invoking utility " << utilityName << endl;

        stringList args(arguments.length() + 1);


        // Convert args to Foam type and add application
        unsigned int argi = 0;

        args[argi++] = utilityName;
        for (unsigned int i = 0; i <arguments.length(); i++)
        {
            args[argi++] = static_cast<const char*>(arguments[i]);
        }

        if (hostName == Foam::hostName())
        {
            // Running on same host.
            pid = ProcessControl::fork(args, logName);

            if (!background && pid > 0)
            {
                // Wait for process to finish.
                ProcessControl::waitpid(pid);
            }
        }
        else
        {
            // Running on different host. No easy way of obtaining pid so
            // use plain system.

            int err = procControl_.system
            (
                procControl_.remoteShellArgs
                (
                    userName(),
                    hostName,
                    args,
                    string(logName),
                    background
                ),
                background ? procControl_.timeOut() : 0
            );

            if (err != 0)
            {
                string msg("Error executing utility ");
                msg += utilityName;

                if (err == ProcessControl::TIMEDOUT)
                {
                    // Timeout will only occur if running in foreground.
                    // (something seriously wrong if spawn does not work in
                    // timeout)
                    msg +=
                        "\nHost "
                        + string(hostName)
                        + " can not be reached (timeout).";

                    if (!background)
                    {
                        throw FoamXSYSError
                        (
                            E_FAIL,
                            msg,
                            hostName,
                            functionName,
                            __FILE__, __LINE__
                        );
                    }
                }
                if (err == -1)
                {
                    msg += "\nCall to system() returned -1";
                }
                else if (err == 127)
                {
                    msg += "\nCall to system() returned 127";
                }
                else
                {
                    msg +=
                        "\nThe utility exited with error code "
                      + name(err);
                }

                throw FoamXError
                (
                    E_FAIL,
                    msg,
                    functionName,
                    __FILE__, __LINE__
                );
            }
        }
    }
    CATCH_ALL(functionName);

    return pid;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::openCase
(
    const FoamXServer::CaseDescriptor& caseDesc
)
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::openCase"
        "(const FoamXServer::CaseDescriptor& caseDesc)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        caseIsInError(caseDesc);

        fileName expandedRootDir = fileName(caseDesc.rootDir).expand();
        fileName caseDir = expandedRootDir/fileName(caseDesc.caseName);
        fileName caseFileName = caseDir/"system/controlDict";
        fileName lockFileName = caseDir/".fxLock";

        if
        (
            orb_.isObjectBound
            (
                caseServerKey(caseDesc.rootDir, caseDesc.caseName)
            )
        )
        {
            throw FoamXError
            (
                E_FAIL,
                "Case " + caseDir
              + " is already registered with nameserver",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Check that the control dictionary exists.
        if (!exists(caseFileName))
        {
            throw FoamXError
            (
                E_FAIL,
                "Case control dictionary not found at '" +  caseFileName + "'",
                functionName,
                __FILE__, __LINE__
            );
        }

        // See if this case is being pre-processed by someone
        // (or something?) else.
        if (exists(lockFileName))
        {
            throw FoamXError
            (
                E_FAIL,
                "Case " + caseDir + " locked",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Everything is pucker so run the case server process.
        log << "Spawning CaseServer process to open the case " 
            << caseDir << endl;

        stringList serverArgs(args_.size() + 3);
        label argI = 0;

        serverArgs[argI++] = "FoamXCaseServer";

        // Pass through caseBrowser args
        for(label i = 1; i < args_.size(); i++)
        {
            serverArgs[argI++] = args_[i];
        }
        serverArgs[argI++] =
            "-case " + expandedRootDir/fileName(caseDesc.caseName);
        serverArgs[argI++] = "-open";
        serverArgs[argI++] = "&";

        if (procControl_.system(serverArgs, procControl_.timeOut()) != 0)
        {
            throw FoamXError
            (
                E_FAIL,
                "Error executing CaseServer process",
                functionName,
                __FILE__, __LINE__
            );
        }
    }
    CATCH_ALL(functionName);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::openCasePost
(
    const FoamXServer::CaseDescriptor& caseDesc,
    const CORBA::Long nProcs
)
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::openCasePost"
        "(const FoamXServer::CaseDescriptor& caseDesc,"
        " const CORBA::Long nProcs)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        fileName expandedRootDir = fileName(caseDesc.rootDir).expand();
        fileName caseDir = expandedRootDir/fileName(caseDesc.caseName);
        fileName caseFileName = caseDir/"system/controlDict";
        fileName lockFileName = caseDir/".fxLock";

        // Check that the control dictionary exists.
        if (!exists(caseFileName))
        {
            throw FoamXError
            (
                E_FAIL,
                "Case control dictionary not found at '" +  caseFileName + "'",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Everything is pucker so run the case server process.
        log << "Spawning CasePostServer process to postprocess the case " 
            << caseDir << endl;

        stringList serverArgs;

        label argI = 0;

        if (nProcs == 1)
        {
            serverArgs.setSize(args_.size() + 3);
        }
        else
        {
            serverArgs.setSize(args_.size() + 5);
            serverArgs[argI++] = "foamJob";
            serverArgs[argI++] = "-ofd";
        }

        serverArgs[argI++] = "FoamXCasePostServer";

        // Pass through caseBrowser args
        for(label i = 1; i < args_.size(); i++)
        {
            serverArgs[argI++] = args_[i];
        }

        serverArgs[argI++] =
            "-case " + expandedRootDir/fileName(caseDesc.caseName);
        serverArgs[argI++] = "-post";
        serverArgs[argI++] = "&";

        Info<< "Spawning:" << serverArgs << endl;

        if (procControl_.system(serverArgs, procControl_.timeOut()) != 0)
        {
            throw FoamXError
            (
                E_FAIL,
                "Error executing CasePostServer process",
                functionName,
                __FILE__, __LINE__
            );
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::newCase
(
    const char* rootDir,
    const char* caseName,
    const char* app
)
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::newCase"
        "(const char* rootDir, const char* caseName, const char* app)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        fileName expandedRootDir = fileName(rootDir).expand();
        fileName caseDir = expandedRootDir/caseName;
        fileName caseFileName = caseDir/"system/controlDict";
        fileName lockFileName = caseDir/".fxLock";

        if (orb_.isObjectBound(caseServerKey(rootDir, caseName)))
        {
            throw FoamXError
            (
                E_FAIL,
                "Case " + caseDir + " is already registered with nameserver",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Check that the control dictionary doesn't exist.
        if (dir(caseDir))
        {
            throw FoamXError
            (
                E_FAIL,
                "Directory already exists at '"
              + caseDir + "'",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Everything is pucker so run the case server process.
        log << "Spawning CaseServer process to create the case " 
            << caseDir << " for application " << app << endl;

        stringList serverArgs(args_.size() + 4);

        label argI = 0;

        serverArgs[argI++] = "FoamXCaseServer";

        // Pass through caseBrowser args
        for(label i = 1; i < args_.size(); i++)
        {
            serverArgs[argI++] = args_[i];
        }

        serverArgs[argI++] = "-case " + expandedRootDir/fileName(caseName);
        serverArgs[argI++] = "-create";
        serverArgs[argI++] = app;
        serverArgs[argI++] = "&";

        if (procControl_.system(serverArgs, procControl_.timeOut()) != 0)
        {
            throw FoamXError
            (
                E_FAIL,
                "Error executing CaseServer process",
                functionName,
                __FILE__, __LINE__
            );
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::importCase
(
    const char* rootDir,
    const char* caseName,
    const char* app
)
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::importCase"
        "(const char* rootDir, const char* caseName, const char* app)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        fileName expandedRootDir = fileName(rootDir).expand();
        fileName caseDir = expandedRootDir/caseName;
        fileName caseFileName = caseDir/"system/controlDict";
        fileName lockFileName = caseDir/".fxLock";

        if
        (
            orb_.isObjectBound
            (
                caseServerKey(rootDir, caseName)
            )
        )
        {
            throw FoamXError
            (
                E_FAIL,
                "Case " + caseDir + " is already registered with nameserver",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Check that the control dictionary exists.
        if (!exists(caseFileName))
        {
            throw FoamXError
            (
                E_FAIL,
                "Case control dictionary not found at '" +  caseFileName + "'",
                functionName,
                __FILE__, __LINE__
            );
        }

        // See if this is being pre-processed by someone (or something?) else.
        if (exists(lockFileName))
        {
            throw FoamXError
            (
                E_FAIL,
                "Case locked",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Everything is pucker so run the case server process.
        log << "Spawning CaseServer process to import the case " 
            << caseDir << " for application " << app << endl;

        stringList serverArgs(args_.size() + 4);

        label argI = 0;

        serverArgs[argI++] = "FoamXCaseServer";

        // Pass through caseBrowser args
        for(label i = 1; i < args_.size(); i++)
        {
            serverArgs[argI++] = args_[i];
        }

        serverArgs[argI++] = "-case " + expandedRootDir/fileName(caseName);
        serverArgs[argI++] = "-import";
        serverArgs[argI++] = app;
        serverArgs[argI++] = "&";

        if (procControl_.system(serverArgs, procControl_.timeOut()) != 0)
        {
            throw FoamXError
            (
                E_FAIL,
                "Error executing CaseServer process",
                functionName,
                __FILE__, __LINE__
            );
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::deleteCase
(
    const FoamXServer::CaseDescriptor& caseDesc
)
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::deleteCase"
        "(const FoamXServer::CaseDescriptor& caseDesc)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        fileName expandedRootDir = fileName(caseDesc.rootDir).expand();
        fileName caseDir = expandedRootDir/fileName(caseDesc.caseName);
        fileName caseFileName = caseDir/"system/controlDict";
        fileName lockFileName = caseDir/".fxLock";

        // Check that the control dictionary exists.
        if (!exists(caseDir))
        {
            throw FoamXError
            (
                E_FAIL,
                "Case control dictionary not found at '" +  caseDir + "'",
                functionName,
                __FILE__, __LINE__
            );
        }

        // See if this is being pre-processed by someone (or something?) else.
        if (exists(lockFileName))
        {
            throw FoamXError
            (
                E_FAIL,
                "Case locked",
                functionName,
                __FILE__, __LINE__
            );
        }

        log << "Deleting case " << caseDir << endl;

        // Delete all the case dictionaries and data files.
        rmDir(caseDir);

        // Remove case from list.
        Foam::HashPtrTable<FoamXServer::CaseDescriptor, string>::iterator
            iter = cases_.find(caseDir);

        if (iter != cases_.end())
        {
            cases_.erase(iter);
        }
        else
        {
            throw FoamXError
            (
                E_FAIL,
                "Cannot remove case " + caseDir + " from case list",
                functionName,
                __FILE__, __LINE__
            );
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::cloneCase
(
    const FoamXServer::CaseDescriptor& caseDesc,
    const char* newRootDir,
    const char* newCaseName,
    const char* newAppClassName,
    const char* timeSel
)
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::cloneCase"
        "(const FoamXServer::CaseDescriptor& caseDesc,"
        " const char* newRootDir, const char* newCaseName,"
        " const char* newAppClassName, const char* timeSel)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        fileName expandedRootDir = fileName(caseDesc.rootDir).expand();
        fileName caseDir = expandedRootDir/fileName(caseDesc.caseName);
        fileName caseFileName = caseDir/"system/controlDict";
        fileName expandedNewRootDir = fileName(newRootDir).expand();
        fileName newCaseDir = expandedNewRootDir/newCaseName;

        // Check that the source control dictionary exists.
        if (!exists(caseFileName))
        {
            throw FoamXError
            (
                E_FAIL,
                "Source case control dictionary not found at '"
              + caseFileName + "'",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Check that the source and destination cases aren't the same.
        if (caseDir == newCaseDir)
        {
            throw FoamXError
            (
                E_FAIL,
                "Invalid destination case name",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Check that the destination case does not yet exist.
        if (exists(newCaseDir))
        {
            throw FoamXError
            (
                E_FAIL,
                "Destination directory '"
              + newCaseDir
              + "' already exists",
                functionName,
                __FILE__, __LINE__
            );
        }

        log << "Copying case " << caseDir << " to "<< newCaseDir << endl;

        if
        (
            !mkDir(newCaseDir)
         || !cp(caseDir/"system", newCaseDir)
         || !cp(caseDir/"constant", newCaseDir)
        )
        {
            throw FoamXError
            (
                E_FAIL,
                "Could not copy system and/or constant directory to '"
              + newCaseDir,
                functionName,
                __FILE__, __LINE__
            );
        }

        instantList times(Time::findTimes(caseDir));

        label firstTimeI = 0;
        label lastTimeI = -1;

        word timeSelMode(timeSel);

        if (timeSelMode == "latestTime")
        {
            firstTimeI = times.size() - 1;
            lastTimeI = firstTimeI;
        }
        else if (timeSelMode == "firstTime")
        {
            firstTimeI = 0;
            lastTimeI = firstTimeI;
        }
        else if (timeSelMode == "noTime")
        {
            firstTimeI = labelMax;
            lastTimeI = labelMin;
        }
        else if (timeSelMode == "allTime")
        {
            firstTimeI = 0;
            lastTimeI = times.size() - 1;
        }
        else
        {
            throw FoamXError
            (
                E_FAIL,
                "Illegal time selection "
              + timeSelMode,
                functionName,
                __FILE__, __LINE__
            );
        }

        for
        (
            label timeI = firstTimeI;
            (timeI <= lastTimeI) && (timeI < times.size());
            timeI++
        )
        {
            log << "Copying time directory "
                << caseDir/times[timeI].name() << endl;

            if (!cp(caseDir/times[timeI].name(), newCaseDir))
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Could not copy time directory '"
                  + caseDir/times[timeI].name()
                  + "'",
                    functionName,
                    __FILE__, __LINE__
                );
            }
        }

        // Create a temporary foam database.
        Time db(Time::controlDictName, newRootDir, newCaseName);

        // Read controlDict, change application and save
        IOdictionary controlDict
        (
            IOobject
            (
                Time::controlDictName,
                db.system(),
                db,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        if
        (
            !controlDict.remove("application")
         && !controlDict.remove("applicationClass")
        )
        {
            throw FoamXIOError
            (
                "Could not find application entry in controlDict",
                controlDict.name(),
                controlDict.startLineNumber(),
                controlDict.endLineNumber(),
                functionName,
                __FILE__, __LINE__
            );
        }

        controlDict.add("application", newAppClassName);

        if (!controlDict.regIOobject::write())
        {
            throw FoamXError
            (
                E_FAIL,
                "Could not write controlDict '"
              + controlDict.path()
              + "'",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Add the new case to the case list
        addCase
        (
            expandedNewRootDir.c_str(),
            newRootDir,
            newCaseName,
            newAppClassName
        );

        // Unlock the new case if it is locked
        unlockCase(expandedNewRootDir.c_str(), newCaseName);

        // TODO : Change the root header entry in all case dictionaries.
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

CORBA::Boolean FoamX::ICaseBrowserImpl::caseLocked
(
    const FoamXServer::CaseDescriptor& caseDesc
)
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::caseLocked"
        "(const FoamXServer::CaseDescriptor& caseDesc)";

    LogEntry log(functionName, __FILE__, __LINE__);

    // See if this is being pre-processed by someone (or something?) else.
    fileName lockFileName =
        fileName(caseDesc.rootDir).expand()/fileName(caseDesc.caseName)
       /".fxLock";

    return exists(lockFileName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::unlockCase
(
    const char* rootDir,
    const char* caseName
)
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::unlockCase"
        "(const FoamXServer::CaseDescriptor& caseDesc)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        fileName caseDir = fileName(rootDir)/caseName;
        fileName lockFileName = caseDir/".fxLock";

        log << "Unlocking case " << caseDir << endl;

        if (exists(lockFileName))
        {
            rm(lockFileName);
        }

        Foam::HashPtrTable<FoamXServer::CaseDescriptor, string>::iterator
            iter = cases_.find(caseDir);

        if (iter != cases_.end())
        {
            iter()->locked = false;
        }
        else
        {
            throw FoamXError
            (
                E_FAIL,
                "Cannot unlock case " + caseDir + ", not in case list",
                functionName,
                __FILE__, __LINE__
            );
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::unlockCaseDescriptor
(
    const FoamXServer::CaseDescriptor& caseDesc
)
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::unlockCase"
        "(const FoamXServer::CaseDescriptor& caseDesc)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        unlockCase(caseDesc.rootDir, caseDesc.caseName);
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::CaseDescriptor* FoamX::ICaseBrowserImpl::caseDescriptor
(
    const fileName& rootDir,
    const fileName& rawRootDir,
    const fileName& caseName
)
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::caseDescriptor"
        "(const fileName& rootDir, const fileName& rawRootDir, "
        "const fileName& caseName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    fileName caseDir = fileName(rootDir)/fileName(caseName);

    try
    {
        log << "Considering case " << caseDir << endl;

        CaseDescriptor* caseDescPtr = new FoamXServer::CaseDescriptor();

        if (caseDescPtr == NULL)
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Failed to create case descriptor for case '" + caseDir,
                functionName,
                __FILE__, __LINE__
            );
        }

        CaseDescriptor& caseDesc = *caseDescPtr;

        // Fill the caseDescriptor
        caseDesc.rootDir = rootDir.c_str();
        caseDesc.rawRootDir = rawRootDir.c_str();
        caseDesc.caseName = caseName.c_str();
        caseDesc.app = "";


        // Get number of processors
        label nProcs = 0;
        while(exists(fileName(caseDir/"processor" + name(nProcs))))
        {
            nProcs++;
        }
        if (nProcs == 0)
        {
            nProcs = 1;
        }
        caseDesc.nProcs = nProcs;


        fileName lockFile = caseDir/".fxLock";

        // See if this case is locked by another user.
        caseDesc.locked = exists(lockFile);

        // Set the error state to false
        caseDesc.error = false;

        // See if the lock file contains an error state
        if (caseDesc.locked)
        {
            // Extract the IOR of the FoamX case server
            IFstream lockFileStream(lockFile);
            word ior(lockFileStream);

            if (ior == "FatalError")
            {
                //reThrow(lockFileStream);
                caseDesc.error = true;
            }
        }

        // See if this is a FoamX managed case
        // (ie, the application entry in the
        // controlDict is defined).
        fileName controlDictFileName = caseDir/"system/controlDict";

        dictionary controlDict((IFstream(controlDictFileName)()));

        if (controlDict.found("application"))
        {
            caseDesc.managed = true;
            word app = controlDict.lookup("application");
            caseDesc.app = app.c_str();
        }
        else if (controlDict.found("applicationClass"))
        {
            caseDesc.managed = true;
            word app = controlDict.lookup("applicationClass");
            caseDesc.app = app.c_str();
        }

        log << "Found case " << caseDir << endl;

        return caseDescPtr;
    }
    CATCH_ALL(functionName);

    return NULL;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::addCase
(
    const char* rootDir,
    const char* rawRootDir,
    const char* caseName,
    const char* appName
)
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::addCase"
        "(const fileName& rootDir, const fileName& rawRootDir, "
        "const fileName& caseName, const fileName& appName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        fileName caseDir = fileName(rootDir)/fileName(caseName);

        log << "Adding case " << caseDir << endl;

        CaseDescriptor* caseDescPtr = new FoamXServer::CaseDescriptor();
        CaseDescriptor& caseDescriptor = *caseDescPtr;

        caseDescriptor.rootDir = rootDir;
        caseDescriptor.rawRootDir = rawRootDir;
        caseDescriptor.caseName = caseName;
        caseDescriptor.app = appName;
        caseDescriptor.managed = true;
        caseDescriptor.locked = true;
        caseDescriptor.error = false;

        if (!cases_.insert(caseDir, caseDescPtr))
        {
            delete caseDescPtr;
            caseDescPtr = NULL;

            throw FoamXError
            (
                E_UNEXPECTED,
                "Failed to add case '" + caseDir
                + "' to case list.  It is probably already in list",
                functionName,
                __FILE__, __LINE__
            );
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::caseOpen
(
    const char* rootDir,
    const char* caseName
)
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::caseOpen"
        "(const fileName& rootDir, const fileName& caseName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        fileName caseDir = fileName(rootDir)/fileName(caseName);

        // Remove case from list.
        Foam::HashPtrTable<FoamXServer::CaseDescriptor, string>::iterator
            iter = cases_.find(caseDir);

        if (iter != cases_.end())
        {
            iter()->locked = true;
            iter()->error = false;
        }
        else
        {
            throw FoamXError
            (
                E_FAIL,
                "Cannot find case " + caseDir + " in case list",
                functionName,
                __FILE__, __LINE__
            );
        }
    }
    CATCH_ALL(functionName);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

CORBA::Boolean FoamX::ICaseBrowserImpl::isCaseInError
(
    const FoamXServer::CaseDescriptor& caseDesc
)
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::isCaseInError"
        "(const FoamXServer::CaseDescriptor& caseDesc)";

    LogEntry log(functionName, __FILE__, __LINE__);

    // See if this is being pre-processed by someone (or something?) else.
    fileName lockFile = 
        fileName(caseDesc.rootDir).expand()/fileName(caseDesc.caseName)
       /".fxLock";

    // See if the lock file contains an error state
    if (exists(lockFile))
    {
        // Extract the IOR of the FoamX case server
        IFstream lockFileStream(lockFile);
        word ior(lockFileStream);

        if (ior == "FatalError")
        {
            reThrow(lockFileStream);
            return true;
        }
    }

    return false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::caseIsInError
(
    const FoamXServer::CaseDescriptor& caseDesc
)
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::caseIsInError"
        "(const FoamXServer::CaseDescriptor& caseDesc)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        fileName caseDir =
            fileName(caseDesc.rootDir)/fileName(caseDesc.caseName);

        // Remove case from list.
        Foam::HashPtrTable<FoamXServer::CaseDescriptor, string>::iterator
            iter = cases_.find(caseDir);

        if (iter != cases_.end())
        {
            iter()->locked = true;
            iter()->error = true;
        }
        else
        {
            throw FoamXError
            (
                E_FAIL,
                "Cannot find case " + caseDir + " in case list",
                functionName,
                __FILE__, __LINE__
            );
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::refreshCaseList()
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::refreshCaseList()";

    LogEntry log(functionName, __FILE__, __LINE__);

    // Collect used roots
    HashTable<fileName> usedRoots;

    for
    (
        Foam::HashPtrTable<FoamXServer::CaseDescriptor, string>::iterator
            iter = cases_.begin();
        iter != cases_.end();
        ++iter
    )
    {
        fileName fName(iter()->rootDir);

        if (!usedRoots.found(fName))
        {
            usedRoots.insert(fName, fName);
        }
    }

    // Clear the existing list
    cases_.clear();

    forAll(usedRoots.toc(), rootI)
    {
        addToCaseList(usedRoots.toc()[rootI].c_str());
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::addToCaseList(const char* dir)
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::refreshCaseList(const char*)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Get case root list.
        StringList_var rootDirs = foamProperties_->rootDirectories();
        StringList_var rawRootDirs = foamProperties_->rawRootDirectories();

        // Loop over all case roots
        unsigned int rootIndex = 0;
        bool rootFound = false;

        for (unsigned int i = 0; i < rootDirs->length(); i++  )
        {
            fileName storedRootDir(rootDirs[i]);

            if (exists(storedRootDir) && (storedRootDir == dir))
            {
                rootFound = true;
                rootIndex = i;
                break;
            }
        }
        if (!rootFound)
        {
            throw FoamXError
            (
                E_FAIL,
                "Cannot find root " + fileName(dir) + " in list of roots",
                functionName,
                __FILE__, __LINE__
            );
        }

        fileName rootDir(rootDirs[rootIndex]);
        fileName rawRootDir(rawRootDirs[rootIndex]);

        // Get a list of directories under this root.
        fileNameList caseDirs
        (
            readDir(rootDir, fileName::DIRECTORY)
        );

        // Search each potential case for a control dictionary.
        for (int nCaseDir = 0; nCaseDir < caseDirs.size(); nCaseDir++)
        {
            try
            {
                fileName controlDictFileName =
                    rootDir/caseDirs[nCaseDir]/"system/controlDict";

                if (exists(controlDictFileName))
                {
                    CaseDescriptor* caseDescPtr = caseDescriptor
                    (
                        rootDir,
                        fileName(rawRootDir),
                        caseDirs[nCaseDir]
                    );

                    if (caseDescPtr)
                    {
                        cases_.insert
                        (
                            rootDir/caseDirs[nCaseDir],
                            caseDescPtr
                        );
                    }
                    else
                    {
                        throw FoamXError
                        (
                            E_UNEXPECTED,
                            "Failed to create caseDescriptor for case '"
                          + rootDir/caseDirs[nCaseDir] + "'",
                            functionName,
                            __FILE__, __LINE__
                        );
                    }
                }
            }
            catch (error& fErr)
            {
                log << fErr << " failed to open case "
                    << rootDir/caseDirs[nCaseDir] << endl;
            }
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::label FoamX::ICaseBrowserImpl::getEntry
(
    const dictionary& dict,
    const word& key,
    const label defaultValue
)
{
    if (dict.found(key))
    {
        return readLabel(dict[key]);
    }
    else
    {
        return defaultValue;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::string FoamX::ICaseBrowserImpl::getEntry
(
    const dictionary& dict,
    const word& key,
    const string& defaultValue
)
{
    if (dict.found(key))
    {
        return dict[key];
    }
    else
    {
        return defaultValue;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::readJobs
(
    JobHashTable& jobsList,
    const fileName& dir,
    const bool throwOnError
)
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::readJobs"
        "(JobHashTable& jobsList, const fileName& dir)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Clear the existing list
        jobsList.clear();

        // Read the names of all the job files
        fileNameList jobFiles = readDir(dir, fileName::FILE);

        forAll (jobFiles, i)
        {
            const fileName& jobFile = jobFiles[i];

            try
            {
                if (size(dir/jobFile) > 0)
                {
                    dictionary jobDict(IFstream(dir/jobFile)());

                    // Required entries

                    JobDescriptor* jobDescPtr = new JobDescriptor();

                    string::size_type dotPos =
                        jobFile.find_last_of('.');

                    jobDescPtr->jobID.hostName = jobFile(0, dotPos).c_str();

                    jobDescPtr->jobID.processID =
                        atoi(jobFile(dotPos + 1, string::npos).c_str());

                    jobDescPtr->startDate =
                        getEntry(jobDict, "startDate", "").c_str();

                    jobDescPtr->startTime = 
                        getEntry(jobDict, "startTime", "").c_str();


                    // Cannot use getEntry since returns string which is
                    // invalid word.
                    word foamUser("");
                    if (jobDict.found("userName"))
                    {
                        foamUser = word(jobDict["userName"]);
                    }
                    jobDescPtr->userName = foamUser.c_str();

                    word foamVersion("0.0");
                    if (jobDict.found("foamVersion"))
                    {
                        foamVersion = token(jobDict["foamVersion"]).wordToken();
                    }
                    jobDescPtr->foamVersion = foamVersion.c_str();

                    word foamCode("");
                    if (jobDict.found("code"))
                    {
                        foamCode = word(jobDict["code"]);
                    }
                    jobDescPtr->code = foamCode.c_str();

                    jobDescPtr->argList = 
                        getEntry(jobDict, "argList", "").c_str();

                    jobDescPtr->currentDir = 
                        getEntry(jobDict, "currentDir", "").c_str();

                    jobDescPtr->pgid =
                        getEntry(jobDict, "PGID", 0);

                    jobDescPtr->ppid = 0;
                        getEntry(jobDict, "PPID", 0);

                    fileName root =
                        getEntry(jobDict, "root", "");

                    if ((root.size() > 0) && (root[0] != '/'))
                    {
                        jobDescPtr->rootDir =
                        (
                            fileName(jobDescPtr->currentDir)
                          + "/" + root
                        ).c_str();
                    }
                    else
                    {
                        jobDescPtr->rootDir =
                            root.c_str();
                    }    

                    jobDescPtr->caseName =
                        getEntry(jobDict, "case", "").c_str();

                    jobDescPtr->nProcs =
                        getEntry(jobDict, "nProcs", 1);

                    if (jobDict.found("slaves"))
                    {
                        wordList slaveJobs(jobDict.lookup("slaves"));
                        jobDescPtr->slaves.length(slaveJobs.size());

                        forAll(slaveJobs, i)
                        {
                            const word& sl = slaveJobs[i];

                            string::size_type dotPos =
                                sl.find_last_of('.');

                            jobDescPtr->slaves[i].hostName =
                                sl(0, dotPos).c_str();
                            jobDescPtr->slaves[i].processID =
                                atoi(sl(dotPos + 1, string::npos).c_str());
                        }
                    }
                    else
                    {
                        jobDescPtr->slaves.length(0);
                    }

                    jobDescPtr->nCountedProcs =
                        getEntry(jobDict, "nCountedProcs", 0);

                    if (jobDict.found("cpuTime"))
                    {
                        jobDescPtr->cpuTime =
                            readScalar(jobDict["cpuTime"]);
                    }
                    else
                    {
                        jobDescPtr->cpuTime = 0.0;
                    }

                    jobDescPtr->endDate = 
                        getEntry(jobDict, "endDate", "").c_str();

                    jobDescPtr->endTime = 
                        getEntry(jobDict, "endTime", "").c_str();

                    jobDescPtr->status = JOB_UNDEFINED;

                    if (jobDict.found("termination"))
                    {
                        word termination(jobDict.lookup("termination"));

                        if (termination == "normal")
                        {
                            jobDescPtr->status = JOB_FINISHED;
                        }
                        else
                        {
                            jobDescPtr->status = JOB_ABORTED;
                        }
                    }

                    jobsList.insert(jobDescPtr->jobID, jobDescPtr);
                }
            }
            catch (...)
            {
                fileName fName(dir/jobFile);


                if (throwOnError)
                {
                    throw FoamXIOError
                    (
                       "Error in job description file",
                        fName,
                        -1, -1,
                        functionName,
                        __FILE__, __LINE__
                    );
                }
                else
                {
                    WarningIn(functionName)
                        << "Removing illegal jobInfo file " << fName
                        << endl;

                    // Try removing jobInfo file (try since we might not own it)
                    rm(fName);
                }
            }
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::refreshJobsLists()
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::refreshJobsLists()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        if (!env("FOAM_JOB_DIR"))
        {
            throw FoamXError
            (
                E_FAIL,
                "Cannot get essential environment variable 'FOAM_JOB_DIR'",
                functionName,
                __FILE__, __LINE__
            );
        }

        fileName licenceDir = Foam::getEnv("FOAM_JOB_DIR");

        // Throw error for illegal job entry files. Don't remove since might
        // upset licensing.
        readJobs(runningJobs_, licenceDir/"runningJobs", true);
        // Remove any illegal job entry files
        readJobs(finishedJobs_, licenceDir/"finishedJobs", false);
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::checkRunningJobs()
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::checkRunningJobs()";

    LogEntry log(functionName, __FILE__, __LINE__);

    fileName licenceDir = Foam::getEnv("FOAM_JOB_DIR");


    // Get hashtable of all hosts running jobs

    HashSet<word> machines;

    try
    {
        if (!env("FOAM_JOB_DIR"))
        {
            throw FoamXError
            (
                E_FAIL,
                "Cannot get essential environment variable 'FOAM_JOB_DIR'",
                functionName,
                __FILE__, __LINE__
            );
        }


        for
        (
            JobHashTable::const_iterator iter = runningJobs_.begin();
            iter != runningJobs_.end();
            ++iter
        )
        {
            machines.insert(word(iter.key().hostName));
        }
    }
    CATCH_ALL(functionName);


    try
    {

        // Run external utility foamProcessInfo to get processes running on
        // hosts. It writes file (see processLogFileName below) which is gives
        // for every process the status.

        stringList args(2);
        args[0] = "foamProcessInfo";


        HashSet<word> unLicensedMachines;
        string unLicensedMsg = "";

        for
        (
            HashSet<word>::const_iterator hostIter = machines.begin();
            hostIter != machines.end();
            ++hostIter
        )
        {
            const word& hostName = hostIter.key();

            log << "Checking host " << hostName << endl;

            try
            {
                if (!hostBrowser_->isHostAlive(hostName.c_str()))
                {
                    log << "Skipping non-accessible host " << hostName << endl;
                    continue;
                }
            }
            catch (FoamXServer::FoamXError& fxErr)
            {
                log << "Non-licensed host " << hostName << endl;

                unLicensedMachines.insert(hostName);

                unLicensedMsg += hostName + " ";

                continue;
            }

            fileName processLogFileName =
                licenceDir/(userName() + '@' + hostName + ".log");

            args[1] = processLogFileName;

            stringList systemArgs = procControl_.remoteShellArgs
            (
                userName(),
                hostName,
                args,
                "",
                false
            );
            int sysStatus = procControl_.system
            (
                systemArgs,
                procControl_.timeOut()
            );

            if (sysStatus == ProcessControl::TIMEDOUT)
            {
                string msg =
                    "Error executing "
                    + procControl_.commandString(systemArgs)
                    + "\nHost " + hostName + " can not be reached (timeout).";

                // Tell Java that machine is non-accessible
                throw FoamXSYSError
                (
                    E_FAIL,
                    msg,
                    hostName,
                    functionName,
                    __FILE__, __LINE__
                );
            }
            else if (sysStatus != 0)
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Error executing "
                        + procControl_.commandString(systemArgs),
                    functionName,
                    __FILE__, __LINE__
                );
            }

            IFstream processLogFile(processLogFileName);

            if (!processLogFile)
            {
                throw FoamXIOError
                (
                    "Error opening file",
                    processLogFileName,
                    -1, -1,
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Convert output from foamProcessInfo into hashtable
            HashTable<word, label, Hash<label> > processTable(processLogFile);
            rm(processLogFileName);

            // Check for each process in runningJobs whether it is
            // in the processLogFile

            for
            (
                JobHashTable::iterator jobIter = runningJobs_.begin();
                jobIter != runningJobs_.end();
                ++jobIter
            )
            {
                if (word(jobIter.key().hostName) == hostName)
                {
                    // See if this process is mentioned in the processTable

                    HashTable<word, label, Hash<label> >::const_iterator
                        processFnd = processTable.find(jobIter.key().processID);

                    if (processFnd != processTable.end())
                    {
                        const word& status = processFnd();

                        if (status == "RUNN")
                        {
                            // Preserve 'ending' status
                            if (jobIter()->status != JOB_STOPPING)
                            {
                                jobIter()->status = JOB_RUNNING;
                            }
                        }
                        else if (status == "SUSP")
                        {
                            jobIter()->status = JOB_SUSPENDED;
                        }
                        else if (status == "OTHR")
                        {
                            jobIter()->status = JOB_UNDEFINED;
                        }
                        else
                        {

                            throw FoamXIOError
                            (
                                "Unknown job status " + status + " in file ",
                                processLogFileName,
                                -1, -1,
                                functionName,
                                __FILE__, __LINE__
                            );
                        }
                    }
                    else
                    {
                        jobIter()->status = JOB_ABORTED;
                    }
                }
            }
        }
        if (unLicensedMachines.size() > 0)
        {
            // Purge runningJobs entries belonging to unlicensed hosts.
            for
            (
                JobHashTable::iterator jobIter = runningJobs_.begin();
                jobIter != runningJobs_.end();
                ++jobIter
            )
            {
                if (unLicensedMachines.found(word(jobIter.key().hostName)))
                {
                    jobIter()->status = JOB_ABORTED;
                }
            }
            purgeRunningJobs();

            string msg =
                "Detected unlicensed hosts\n"
              + unLicensedMsg
              + "\nPurged jobs from runningJobs/.";

            throw FoamXError
            (
                E_FAIL,
                msg,
                functionName,
                __FILE__, __LINE__
            );
        }

    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::purgeRunningJobs()
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::purgeRunningJobs()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        if (!env("FOAM_JOB_DIR"))
        {
            throw FoamXError
            (
                E_FAIL,
                "Cannot get essential environment variable 'FOAM_JOB_DIR'",
                functionName,
                __FILE__, __LINE__
            );
        }

        fileName licenceDir = Foam::getEnv("FOAM_JOB_DIR");
        fileName runningJobsDir = licenceDir/"runningJobs";
        fileName finishedJobsDir = licenceDir/"finishedJobs";

        // 1. Mark jobs to be removed. Cannot remove from runningJobs here since
        //    we're looping over it.
        DynamicList<JobID> toBeRemoved;

        for
        (
            JobHashTable::const_iterator iter = runningJobs_.begin();
            iter != runningJobs_.end();
            ++iter
        )
        {
            if (iter()->status == JOB_ABORTED || iter()->status == JOB_FINISHED)
            {
                toBeRemoved.append(iter.key());
            }
        }

        // 2. Actually remove
        forAll(toBeRemoved, i)
        {
            const JobID& jobID = toBeRemoved[i];

            // Remove any identically named entries in finishedJobs (if any)
            JobHashTable::iterator finishFnd = finishedJobs_.find(jobID);

            if (finishFnd != finishedJobs_.end())
            {
                finishedJobs_.erase(finishFnd);
            }

            // Transfer jobId from runningJobs to finishedJobs
            JobHashTable::iterator runFnd = runningJobs_.find(jobID);

            if (runFnd != runningJobs_.end())
            {
                finishedJobs_.insert(jobID, runningJobs_.remove(runFnd));
            }

            fileName jobFileName =
                word(jobID.hostName)
              + '.'
              + name(jobID.processID);

            if (exists(runningJobsDir/jobFileName))
            {
                mv(runningJobsDir/jobFileName, finishedJobsDir);
            }
        }

        // 3. Remove zero sized files
        fileNameList jobFiles = readDir(runningJobsDir, fileName::FILE);

        forAll (jobFiles, i)
        {
            fileName jobFileName = runningJobsDir/jobFiles[i];

            if (size(jobFileName) == 0)
            {
                rm(jobFileName);
            }
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::purgeFinishedJob(const FoamXServer::JobID& jobID)
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::purgeFinishedJob(const FoamXServer::JobID&)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        if (!env("FOAM_JOB_DIR"))
        {
            throw FoamXError
            (
                E_FAIL,
                "Cannot get essential environment variable 'FOAM_JOB_DIR'",
                functionName,
                __FILE__, __LINE__
            );
        }

        fileName finishedJobsDir = Foam::getEnv("FOAM_JOB_DIR")/"finishedJobs";

        fileName jobFileName = 
            finishedJobsDir
           /(word(jobID.hostName) + '.' + name(jobID.processID));

        rm(jobFileName);

        JobHashTable::iterator finishFnd = finishedJobs_.find(jobID);

        if (finishFnd != finishedJobs_.end())
        {
            finishedJobs_.erase(finishFnd);
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::purgeFinishedJobs(CORBA::Long nDays)
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::purgeFinishedJobs(CORBA::Long nDays)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        if (!env("FOAM_JOB_DIR"))
        {
            throw FoamXError
            (
                E_FAIL,
                "Cannot get essential environment variable 'FOAM_JOB_DIR'",
                functionName,
                __FILE__, __LINE__
            );
        }

        fileName finishedJobsDir = Foam::getEnv("FOAM_JOB_DIR")/"finishedJobs";
        fileName runningJobsDir = Foam::getEnv("FOAM_JOB_DIR")/"runningJobs";

        // Current time to compare against file modification time.
        time_t now = Foam::clock::getTime();

        // 1. Mark jobs to be removed. Cannot remove from finishedJobs here
        //    since we're looping over it.
        DynamicList<JobID> toBeRemoved;

        for
        (
            JobHashTable::const_iterator iter = finishedJobs_.begin();
            iter != finishedJobs_.end();
            ++iter
        )
        {
            fileName jobName =
                word(iter.key().hostName)
              + '.'
              + name(iter.key().processID);


            fileName jobFileName = finishedJobsDir/jobName;

            if (exists(jobFileName))
            {
                if (size(jobFileName) > 0)
                {
                    if ((now - lastModified(jobFileName))/86400 > nDays)
                    {
                        // File older than nDays days
                        rm(jobFileName);
                        toBeRemoved.append(iter.key());
                    }
                    else
                    {
                        //keep file
                    }
                }
                else
                {
                    // Zero sized file. Keep if has runningJobs equivalence
                    if (!exists(runningJobsDir/jobName))
                    {
                        rm(jobFileName);
                        toBeRemoved.append(iter.key());
                    }
                }
            }
            else
            {
                // Does not exist in finishedJobs anymore
                toBeRemoved.append(iter.key());
            }
        }


        // 2. Actually remove
        forAll(toBeRemoved, i)
        {
            JobHashTable::iterator finishFnd =
                finishedJobs_.find(toBeRemoved[i]);

            if (finishFnd != finishedJobs_.end())
            {
                finishedJobs_.erase(finishFnd);
            }
        }

        // 3. Zero size files not in finishedJobs
        fileNameList jobFiles = readDir(finishedJobsDir, fileName::FILE);

        forAll (jobFiles, i)
        {
            fileName jobFileName = finishedJobsDir/jobFiles[i];

            if (size(jobFileName) == 0)
            {
                // Zero sized file. Keep if has runningJobs equivalence
                if (!exists(runningJobsDir/jobFiles[i]))
                {
                    rm(jobFileName);
                }
            }
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::kill(const JobID& jobID)
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::kill(const JobID&)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        int retVal = procControl_.kill(word(jobID.hostName), jobID.processID);

        if (retVal == ProcessControl::TIMEDOUT)
        {
            string msg =
                "Error killing job on " + word(jobID.hostName)
                + " with PID " + name(jobID.processID)
                + "\nHost " + word(jobID.hostName)
                + " can not be reached (timeout).";

            throw FoamXSYSError
            (
                E_FAIL,
                msg,
                word(jobID.hostName),
                functionName,
                __FILE__, __LINE__
            );
        }
        else if (retVal == -1)
        {
            throw FoamXError
            (
                E_FAIL,
                "Error killing job on " + word(jobID.hostName)
              + " with PID " + name(jobID.processID),
                functionName,
                __FILE__, __LINE__
            );
        }
        else
        {
            JobHashTable::iterator iter = runningJobs_.find(jobID);

            if (iter != runningJobs_.end())
            {
                // Note: make copy of key. Needed because key gets destroyed
                // in remove below.
                JobID job(iter.key());

                iter()->status = JOB_ABORTED;

                // Remove from finishedJobs if already present
                JobHashTable::iterator finishFnd = finishedJobs_.find(jobID);

                if (finishFnd != finishedJobs_.end())
                {
                    finishedJobs_.erase(finishFnd);
                }

                // Transfer from runningJobs to finishedJobs.
                finishedJobs_.insert(job, runningJobs_.remove(iter));
            }
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::suspend(const JobID& jobID)
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::suspend(const JobID&)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        int retVal = procControl_.suspend
        (
            word(jobID.hostName),
            jobID.processID
        );
        if (retVal  == ProcessControl::TIMEDOUT)
        {
            string msg =
                "Error killing job on " + word(jobID.hostName)
                + " with PID " + name(jobID.processID)
                + "\nHost " + word(jobID.hostName)
                + " can not be reached (timeout).";

            throw FoamXSYSError
            (
                E_FAIL,
                msg,
                word(jobID.hostName),
                functionName,
                __FILE__, __LINE__
            );
        }
        else if (retVal == -1)
        {
            throw FoamXError
            (
                E_FAIL,
                "Error suspending job on " + word(jobID.hostName)
              + " with PID " + name(jobID.processID),
                functionName,
                __FILE__, __LINE__
            );
        }
        else
        {
            JobHashTable::iterator iter = runningJobs_.find(jobID);

            if (iter != runningJobs_.end())
            {
                iter()->status = JOB_SUSPENDED;
            }
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::cont(const JobID& jobID)
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::cont(const JobID&)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        int retVal = procControl_.cont(word(jobID.hostName), jobID.processID);

        if (retVal == ProcessControl::TIMEDOUT)
        {
            string msg =
                "Error killing job on " + word(jobID.hostName)
                + " with PID " + name(jobID.processID)
                + "\nHost " + word(jobID.hostName)
                + " can not be reached (timeout).";

            throw FoamXSYSError
            (    
                E_FAIL,
                msg,
                word(jobID.hostName),
                functionName,
                __FILE__, __LINE__
            );
        }
        if (retVal == -1)
        {
            throw FoamXError
            (    
                E_FAIL,
                "Error continuing job on " + word(jobID.hostName)
              + " with PID " + name(jobID.processID),
                functionName,
                __FILE__, __LINE__
            );
        }
        else
        {
            JobHashTable::iterator iter = runningJobs_.find(jobID);

            if (iter != runningJobs_.end())
            {
                iter()->status = JOB_RUNNING;
            }
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::end
(
    const JobID& jobID,
    const char* rootDir,
    const char* caseName,
    const CORBA::Boolean now
)
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::end"
        "(const JobID&, const char*, const char*, CORBA::Boolean)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        int argsLen = 3;

        if (now)
        {
            argsLen = 4;
        }

        stringList args(argsLen);

        int argi = 0;

        args[argi++] = "foamEndJob";
        if (now)
        {
            args[argi++] = "-n";
        }
        args[argi++] = "-case " + fileName(rootDir)/fileName(caseName);
        args[argi++] = name(jobID.processID);

        stringList systemArgs = procControl_.remoteShellArgs
        (
            userName(),
            word(jobID.hostName),
            args,
            "",
            true
        );

        int sysStatus = procControl_.system
        (
            systemArgs,
            procControl_.timeOut()
        );

        if (sysStatus == ProcessControl::TIMEDOUT)
        {
            string msg =
                "Error executing "
                + procControl_.commandString(systemArgs)
                + "\nHost " + word(jobID.hostName)
                + " can not be reached (timeout).";

            throw FoamXSYSError
            (
                E_FAIL,
                msg,
                word(jobID.hostName),
                functionName,
                __FILE__, __LINE__
            );
        }
        else if (sysStatus != 0)
        {
            throw FoamXError
            (
                E_FAIL,
                "Error executing "
                    + procControl_.commandString(systemArgs),
                functionName,
                __FILE__, __LINE__
            );
        }
        else
        {
            JobHashTable::iterator iter = runningJobs_.find(jobID);

            if (iter != runningJobs_.end())
            {
                iter()->status = JOB_STOPPING;
            }
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::setStatus
(
    const FoamXServer::JobID& jobID,
    FoamXServer::JobStatus jobStatus
)
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::setStatus"
        "(const FoamXServer::JobID&, FoamXServer::JobStatus)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        JobHashTable::iterator iter = runningJobs_.find(jobID);

        if (iter != runningJobs_.end())
        {
            iter()->status = jobStatus;
        }
        else
        {
            throw FoamXError
            (    
                E_FAIL,
                "Error setting status for job on " + word(jobID.hostName)
              + " with PID " + name(jobID.processID)
              + ", job not in runningJobs list.",
                functionName,
                __FILE__, __LINE__
            );
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::close()
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::close()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Release Foam System Properties object.
        if (foamProperties_ != NULL)
        {
            foamProperties_->_remove_ref();
        }

        // Unregister object from name server.
        NameServer fxNameServer(orb_.orbPtr_);
        fxNameServer.unbindObject(objectName_);

        fxNameServer.removeContext(hostContext_);

        // Explicit disconnect (before orb shutdown)
        fxNameServer.disconnect();

        // Shutdown th'Orb.
        // Removed because it screws-up Java 1.4 CORBA connection
        orb_.orbPtr_->shutdown(false);       // Do not wait for completion.
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::validate()
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::validate()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::ICaseBrowserImpl::save()
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::save()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Get the properties object to save the user properties.
        foamProperties_->saveUserProperties();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::fileName FoamX::ICaseBrowserImpl::caseServerKey
(
    const char* rootDir,
    const char* caseName
)
{
    return hostName()/userName()/fileName(rootDir)/caseName;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

CORBA::Boolean FoamX::ICaseBrowserImpl::getCaseServerReference
(
    const char* rootDir,
    const char* caseName,
    FoamXServer::CaseServer::ICaseServer_out caseObj
)
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::getCaseServerReference"
        "(const char* rootDir, const char* caseName, "
        "FoamXServer::CaseServer::ICaseServer_out caseObj)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        fileName lockFile = fileName(rootDir)/fileName(caseName)/".fxLock";

        if (exists(lockFile))
        {
            // Extract the IOR of the FoamX case server process.
            IFstream lockFileStream(lockFile);
            word ior(lockFileStream);

            if (ior == "FatalError")
            {
                reThrow(lockFileStream);
            }

            log << "IOR : " << ior << endl;

            CORBA::Object_var obj = orb_.orbPtr_->string_to_object(ior.c_str());
            ICaseServer_var ref = ICaseServer::_narrow(obj);

            if (CORBA::is_nil(ref))
            {
                throw FoamXError
                (
                    E_FAIL,
                    "Can't get narrow reference to type ICaseServer "
                    "(or it was nil)",
                    functionName,
                    __FILE__, __LINE__
                );
            }

            // Return reference to the case server object.
            caseObj = ICaseServer::_duplicate(ref);

            return true;
        }
        else
        {
            return false;
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

CORBA::Boolean FoamX::ICaseBrowserImpl::getCasePostServerReference
(
    const char* rootDir,
    const char* caseName,
    CORBA::Long nProcs,
    FoamXServer::CasePostServer::ICasePostServer_out caseObj
)
{
    static const char* functionName =
        "FoamX::ICaseBrowserImpl::getCasePostServerReference"
        "(const char* rootDir, const char* caseName, CORBA::Long nProcs, "
        "FoamXServer::CasePostServer::ICasePostServer_out caseObj)";

    LogEntry log(functionName, __FILE__, __LINE__);

    if (nProcs == 0)
    {
        nProcs = 1;
    }
    fileNameList serverKeys(nProcs);

    serverKeys[0] =
        "FoamXCasePostServer"/userName()
        /string::validate<word>(rootDir)/caseName;

    for(label procIndex = 1; procIndex < nProcs; procIndex++)
    {
        serverKeys[procIndex] =
            "FoamXCasePostServer"/userName()/string::validate<word>(rootDir)
          / word(caseName/("processor" + name(procIndex)));
    }

    try
    {
        // Connect to name server.
        NameServer fxNameServer(orb_.orbPtr_);

        string procLabel = '[' + word(name(Pstream::myProcNo())) + "]-";

        // See if all processors bound
        forAll(serverKeys, keyI)
        {
            if (!fxNameServer.isObjectBound(serverKeys[keyI]))
            {
                Pout<< serverKeys[keyI]
                    << " not yet resolved." << endl;

                return false;
            }
        }

        // Resolve and return master
        caseObj =
        fxNameServer.resolve<FoamXServer::CasePostServer::ICasePostServer>
        (
            serverKeys[0]
        );

        return true;
    }
    CATCH_ALL(functionName);
}


// ************************************************************************* //
