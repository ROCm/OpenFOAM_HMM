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
#include "argList.H"
#include "OSspecific.H"

// FoamX header files.
#include "FoamXErrors.H"
#include "ICaseServerImpl.H"
#include "LogManager.H"
#include "LogEntry.H"
#include "Orb.H"

// Declare all the namespaces used by FoamX
#include "FoamXNameSpaces.H"

using namespace FoamX;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char **argv)
{
    static const char* functionName =
        "FoamXCaseServer::main(int argc, char **argv)";

    fileName lockFileName;

    try
    {
        FatalError.throwExceptions();
        FatalIOError.throwExceptions();

        // Construct the ORB, remove ORB arguments
        Orb orb(argc, argv);

        argList::validOptions.insert("open", "");
        argList::validOptions.insert("create", "app");
        argList::validOptions.insert("import", "app");
        //argList::validOptions.insert("test", "");
        argList args(argc, argv);

        lockFileName = args.path()/".fxLock";

        // Initialise the log manager object. 
        // Write the log into the case directory.
        LogManager logMgr(args.path()/"CaseServerLog.xml");

        // Register the root function call.
        LogEntry log(functionName, __FILE__, __LINE__);

        word mode;
        word app;
        if (args.options().found("open"))
        {
            mode = "open";
        }
        else if (args.options().found("create"))
        {
            mode = "create";
            app = args.options()["create"];
        }
        else if (args.options().found("import"))
        {
            mode = "import";
            app = args.options()["import"];
        }
        else
        {
            FatalErrorIn(functionName)
                << "Essential option not given" << endl;
            args.printUsage();
            FatalError.abort();
        }

        // Allocate the ICaseServerImpl object on the heap.
        // This is a reference counted object and will be deleted by
        // the POA when it is no longer needed.
        ICaseServerImpl* caseServerPtr = new ICaseServerImpl
        (
            orb,
            args.rootPath(),
            args.caseName(),
            mode,
            app
        );

        if (caseServerPtr == NULL)
        {
            FatalErrorIn(functionName)
                << "Failed to instantiate ICaseServerImpl object."
                << exit(FatalError);
        }

        log << "Created ICaseServerImpl object." << endl;

        // Get the object IOR and write to lock file.
        string ior = orb.ior(caseServerPtr->_this());

        log << "IOR : " << ior << "." << endl;

        {
            OFstream ofLock(lockFileName);

            if (ofLock.good())
            {
                ofLock << ior.c_str();
            }
            else
            {
                throw FoamX::FoamXError
                (
                    E_FAIL,
                    "Failed to open lock file " + lockFileName,
                    functionName,
                    __FILE__, __LINE__
                );
            }
        }

        if (args.options().found("test"))
        {
            caseServerPtr->readMesh();
            caseServerPtr->save();
        }
        else
        {
            log << "CaseServer running....." << endl;
            Info << "CaseServer running....." << endl;

            // Go!
            orb.run();

            // Clean up all the resources.
            log << "Shutting down Case Server." << endl;
        }
    }
    catch (Foam::IOerror& fIOErr)
    {
        Serr<< "Caught Foam IOerror in " << functionName << " : "
            << nl << fIOErr << endl;

        writeException(fIOErr, lockFileName);
    }
    catch (Foam::error& fErr)
    {
        Serr<< "Caught Foam error in " << functionName << " : "
            << nl << fErr << endl;

        writeException(fErr, lockFileName);
    }
    catch (FoamXServer::FoamXIOError& fIOErr)
    {
        Serr<< "Caught FoamXIOError exception in " << functionName << " : "
            << nl << fIOErr << endl;

        writeException(fIOErr, lockFileName);
    }
    catch (FoamXServer::FoamXError& ex)
    {
        Serr<< "Caught FoamXError exception in " << functionName << " : "
            << nl << ex << endl;

        writeException(ex, lockFileName);
    }
    catch (CORBA::COMM_FAILURE& ex)
    {
        Serr<< "Caught CORBA::COMM_FAILURE exception in "
            << functionName << endl;

        writeException(ex, lockFileName);
    }
    catch (CORBA::SystemException& ex)
    {
        Serr<< "Caught CORBA::SystemException exception in "
            << functionName << endl;

        writeException(ex, lockFileName);
    }
    catch (CORBA::Exception& ex)
    {
        Serr<< "Caught CORBA::Exception exception in "
            << functionName << endl;

        writeException(ex, lockFileName);
    }
    catch (...)
    {
        Serr<< "Caught system exception in "
            << functionName << endl;

        writeException(systemError(), lockFileName);
    }

    Info<< "Finishing " << functionName << endl;

    return 0;
}


// ************************************************************************* //
