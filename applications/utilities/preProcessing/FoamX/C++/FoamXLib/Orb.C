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

// FoamX header files.
#include "FoamX.H"
#include "FoamXErrors.H"
#include "Orb.H"
#include "NameServer.H"
#include "LogEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::Orb::Orb(int& argc, char** argv)
:
    // Initialise the ORB, passing in the command line arguments.
    orbPtr_(CORBA::ORB_init(argc, argv))
{
    // Obtain a reference to the root POA.
    CORBA::Object_var obj = orbPtr_->resolve_initial_references("RootPOA");
    PortableServer::POA_var poaPtr_ = PortableServer::POA::_narrow(obj);

    // Obtain POAManager reference and tell the POA to start accepting
    // requests on its objects.
    PortableServer::POAManager_var pman = poaPtr_->the_POAManager();
    pman->activate();
}


FoamX::Orb::~Orb()
{
    //poaPtr_->destroy(true, false);
    orbPtr_->destroy();
}


void FoamX::Orb::run()
{
    orbPtr_->run();
}


Foam::string FoamX::Orb::ior(CORBA::Object_ptr oPtr)
{
    return orbPtr_->object_to_string(oPtr);
}


bool FoamX::Orb::isObjectBound(const Foam::fileName& objectName)
{
    static const char* functionName =
        "FoamX::Orb::isObjectBound(const Foam::fileName& objectName)";

    Foam::LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Connect to name server.
        NameServer fxNameServer(orbPtr_);
        return fxNameServer.isObjectBound(objectName);
    }
    CATCH_ALL(functionName);
}


// ************************************************************************* //
