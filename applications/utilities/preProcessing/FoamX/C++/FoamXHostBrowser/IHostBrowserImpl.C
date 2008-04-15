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
#include "dictionary.H"
#include "wordList.H"
#include "stringList.H"
#include "fileNameList.H"
#include "IFstream.H"
#include "OSspecific.H"

// FoamX header files.
#include "FoamXErrors.H"
#include "IHostBrowserImpl.H"
#include "LogEntry.H"
#include "NameServer.H"
#include "Orb.H"
#include "Paths.H"

// Declare all the namespaces used by FoamX.
#include "FoamXNameSpaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IHostBrowserImpl::IHostBrowserImpl
(
    Orb& orb,
    const stringList& args
)
:
    orb_(orb),
    args_(args),
    objectName_("FoamXHostBrowser"),
    procControl_(dotFoam("apps/FoamX/FoamX.cfg"))
{
    static const char* functionName =
        "FoamX::IHostBrowserImpl::IHostBrowserImpl(Orb& orb)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        wordList hosts(1, hostName());

        // Open the user's Foam control dictionary.
        fileName controlDictFileName = dotFoam("controlDict");
        if (exists(controlDictFileName))
        {
            dictionary controlDict((IFstream(controlDictFileName)()));

            // Read hosts list
            if (controlDict.found("hosts"))
            {
                wordList controlDictHosts(controlDict.lookup("hosts"));
                if (controlDictHosts.size())
                {
                    hosts = controlDictHosts;
                }
            }
        }

        forAll(hosts, i)
        {
            FoamXServer::HostDescriptor* hostDescPtr =
                new FoamXServer::HostDescriptor();

            hostDescPtr->name = hosts[i].c_str();
            hostDescPtr->alive = true;

            hosts_.insert(hosts[i], hostDescPtr); 
        }

        // Register this object with the name server.
        // Will throw an exception if any error occurs.
        NameServer fxNameServer(orb_.orbPtr_);
        FoamXServer::HostBrowser::IHostBrowser_var ref = _this();

        // Rebind if necessary.
        fxNameServer.bindObject(objectName_, ref, true);

        // Decrement the reference count of the object implementation, so
        // that it will be properly cleaned up when the POA has determined
        // that it is no longer needed.
        _remove_ref();
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::IHostBrowserImpl::~IHostBrowserImpl()
{
    static const char* functionName =
        "FoamX::IHostBrowserImpl::~IHostBrowserImpl()";

    LogEntry log(functionName, __FILE__, __LINE__);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamXServer::HostDescriptorList* FoamX::IHostBrowserImpl::hosts()
{
    static const char* functionName =
        "FoamXServer::StringList* FoamX::IHostBrowserImpl::hosts()";

    LogEntry log(functionName, __FILE__, __LINE__);


    HostDescriptorList* hostListPtr = new HostDescriptorList();
    hostListPtr->length(hosts_.size());

    label i = 0;
    for
    (
        Foam::HashPtrTable<FoamXServer::HostDescriptor, string>::iterator
            iter = hosts_.begin();
        iter != hosts_.end();
        ++iter
    )
    {
        (*hostListPtr)[i++] = *iter();
    }

    return hostListPtr;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IHostBrowserImpl::refreshHostList()
{
    static const char* functionName =
        "void FoamX::IHostBrowserImpl::refreshHostList()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        for
        (
            Foam::HashPtrTable<FoamXServer::HostDescriptor, string>::
                iterator iter = hosts_.begin();
            iter != hosts_.end();
            ++iter
        )
        {
            try
            {
                iter()->alive =
                    ping(word(iter()->name), procControl_.timeOut());
            }
            catch (error& fErr)
            {
                /* even ping gives problems so rsh/ssh will certainly fail */
                iter()->alive = false;
            }
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

CORBA::Boolean FoamX::IHostBrowserImpl::isHostAlive(const char* hostName)
{
    static const char* functionName =
        "FoamX::IHostBrowserImpl::isHostAlive(const char* hostName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        Foam::HashPtrTable<FoamXServer::HostDescriptor, string>::iterator
            iter = hosts_.find(hostName);

        if (iter != hosts_.end())
        {
            return iter()->alive;
        }
        else
        {
            // Note: don't throw FoamXSYSError since whole host is unknown.
            throw FoamXError
            (
                E_FAIL,
                "Cannot find host '" + word(hostName)
              + "' in list of licenced hosts",
                functionName,
                __FILE__, __LINE__
            );
        }

    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IHostBrowserImpl::hostIsAlive(const char* hostName)
{
    static const char* functionName =
        "FoamX::IHostBrowserImpl::hostIsAlive(const char* hostName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        Foam::HashPtrTable<FoamXServer::HostDescriptor, string>::iterator
            iter = hosts_.find(hostName);

        if (iter != hosts_.end())
        {
            iter()->alive = true;
        }
        else
        {
            throw FoamXError
            (
                E_FAIL,
                "Cannot find host '" + word(hostName)
              + "' in list of licenced hosts",
                functionName,
                __FILE__, __LINE__
            );
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IHostBrowserImpl::hostIsDead(const char* hostName)
{
    static const char* functionName =
        "FoamX::IHostBrowserImpl::hostIsDead(const char* hostName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        Foam::HashPtrTable<FoamXServer::HostDescriptor, string>::iterator
            iter = hosts_.find(hostName);

        if (iter != hosts_.end())
        {
            iter()->alive = false;
        }
        else
        {
            throw FoamXError
            (
                E_FAIL,
                "Cannot find host '" + word(hostName)
              + "' in list of licenced hosts",
                functionName,
                __FILE__, __LINE__
            );
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IHostBrowserImpl::openCaseBrowser
(
    const char* hostName
)
{
    static const char* functionName =
        "FoamX::IHostBrowserImpl::openCaseBrowser"
        "(const char* hostName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Case browser will register itself with the name server under
        // a "hostName/"FoamXCaseBrowser"" key.
        fileName caseBrowserObjectName = fileName(hostName)/"FoamXCaseBrowser";


        if (orb_.isObjectBound(caseBrowserObjectName))
        {
            throw FoamXError
            (
                E_FAIL,
                string("CaseBrowser already running on ") + hostName
              + " according to nameserver",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Start a case browser process on the specified host. Pass through
        // arguments with which HostBrowser was called.
        stringList browserArgs(args_);

        browserArgs[0] = "FoamXCaseBrowser";

        stringList argList
        (
            procControl_.remoteShellArgs
            (
                userName(),
                hostName,
                browserArgs,
                "",         // no log file
                true        // run in background
            )
        );

        if (procControl_.system(argList) != 0)
        {
            hostIsDead(hostName);

            throw FoamXSYSError
            (
                E_FAIL,
                "Server '" + caseBrowserObjectName + "' timed-out.",
                hostName,
                functionName,
                __FILE__, __LINE__
            );
        }
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

CORBA::Boolean FoamX::IHostBrowserImpl::getCaseBrowserReference
(
    const char *hostName,
    FoamXServer::CaseBrowser::ICaseBrowser_out browserObj
)
{
    static const char* functionName =
        "FoamX::IHostBrowserImpl::getCaseBrowserReference"
        "(const char *hostName, "
        "FoamXServer::CaseBrowser::ICaseBrowser_out browserObj)";

    LogEntry log(functionName, __FILE__, __LINE__);

    fileName browserKey(fileName(hostName)/"FoamXCaseBrowser");

    try
    {
        // Connect to name server.
        NameServer fxNameServer(orb_.orbPtr_);

        if (fxNameServer.isObjectBound(browserKey))
        {
            // Return server reference.
            browserObj = 
            fxNameServer.resolve<FoamXServer::CaseBrowser::ICaseBrowser>
            (browserKey);

            return true;
        }
        return false;
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IHostBrowserImpl::validate()
{
    static const char* functionName =
        "FoamX::IHostBrowserImpl::validate()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IHostBrowserImpl::save()
{
    static const char* functionName =
        "FoamX::IHostBrowserImpl::save()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
    }
    CATCH_ALL(functionName);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::IHostBrowserImpl::close()
{
    static const char* functionName =
        "FoamX::IHostBrowserImpl::close()";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Unregister object from the naming service.
        NameServer fxNameServer(orb_.orbPtr_);
        fxNameServer.unbindObject(objectName_);
        fxNameServer.disconnect();

        // Shutdown th'Orb.
        // Remove for Java 1.4 since it messes up CORBA connectionx
        orb_.orbPtr_->shutdown(false);
    }
    CATCH_ALL(functionName);
}


// ************************************************************************* //
