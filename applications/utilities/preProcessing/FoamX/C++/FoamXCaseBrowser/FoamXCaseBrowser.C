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
#include "argList.H"
#include "OSspecific.H"

// FoamX header files.
#include "FoamXErrors.H"
#include "ICaseBrowserImpl.H"
#include "LogManager.H"
#include "LogEntry.H"
#include "Paths.H"
#include "Orb.H"

// Declare all the namespaces used by FoamX
#include "FoamXNameSpaces.H"

using namespace FoamX;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char **argv)
{
    static const char* functionName =
        "FoamXCaseBrowser::main(int argc, char **argv)";

    try
    {
        FatalError.throwExceptions();
        FatalIOError.throwExceptions();

        // Make copy of raw args (since Orb modifies them)
        stringList args(argc);
        forAll(args, argI)
        {
            args[argI] = argv[argI];
        }

        // Initialise the log manager object.
        LogManager logManager(Paths::tmp/"CaseBrowserLog.xml");

        // Register the root function call.
        LogEntry log(functionName, __FILE__, __LINE__);

        // Construct the ORB.
        FoamX::Orb orb(argc, argv);

        // Construct argList so job gets registered.
        argList::validArgs.clear();
        argList dummyArgs(argc, argv);

        // Allocate the ICaseBrowserImpl object on the heap.
        // This is a reference counted object and will be deleted by
        // the POA when it is no longer needed.
        ICaseBrowserImpl* caseBrowserPtr = new ICaseBrowserImpl(orb, args);

        if (caseBrowserPtr == NULL)
        {
            FatalErrorIn(functionName)
                << "Failed to instantiate ICaseBrowserImpl object."
                << exit(FatalError);
        }

        log << "Created ICaseBrowserImpl object" << endl;

        // Emit the IOR.
        log << "IOR : " << orb.ior(caseBrowserPtr->_this()) << endl;

        Info << "CaseBrowser running....." << endl;

        // Go! Go! Go!
        orb.run();
    }
    catch (Foam::IOerror& fIOErr)
    {
        Serr<< "Caught Foam IOerror in " << functionName << " : "
            << nl << fIOErr << endl;
    }
    catch (Foam::error& fErr)
    {
        Serr<< "Caught Foam error in " << functionName << " : "
            << nl << fErr << endl;
    }
    catch (FoamXServer::FoamXIOError& fIOErr)
    {
        Serr<< "Caught FoamXIOError exception in " << functionName << " : "
            << nl << fIOErr << endl;
    }
    catch (FoamXServer::FoamXError& ex)
    {
        Serr<< "Caught FoamXError exception in " << functionName << " : "
            << nl << ex << endl;
    }
    catch (CORBA::COMM_FAILURE& ex)
    {
        Serr<< "Caught CORBA::COMM_FAILURE exception in "
            << functionName << endl;
    }
    catch (CORBA::SystemException& ex)
    {
        Serr<< "Caught CORBA::SystemException exception in "
            << functionName << endl;
    }
    catch (CORBA::Exception& ex)
    {
        Serr<< "Caught CORBA::Exception exception in "
            << functionName << endl;
    }
    catch (...)
    {
        Serr<< "Caught system exception in "
            << functionName << endl;
    }

    Info<< "Finishing " << functionName << endl;

    return 0;
}


// ************************************************************************* //
