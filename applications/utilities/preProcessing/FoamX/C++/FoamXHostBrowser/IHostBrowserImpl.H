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
    IHostBrowserImpl.C

\*---------------------------------------------------------------------------*/

#ifndef IHostBrowserImpl_H
#define IHostBrowserImpl_H

// Foam header files.
#include "HashPtrTable.H"
#include "stringList.H"

// FoamX header files.
#include "FoamX.H"
#include "FoamXWordList.H"
#include "ProcessControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace FoamX
{

// Forward declaration of classes.
class Orb;

/*---------------------------------------------------------------------------*\
                    Class IHostBrowserImpl Declaration
\*---------------------------------------------------------------------------*/

class IHostBrowserImpl
:
    public POA_FoamXServer::HostBrowser::IHostBrowser,
    public PortableServer::RefCountServantBase
{
    // Private data

        //- Reference to the ORB
        Orb& orb_;

        //- Command line arguments
        Foam::stringList args_;

        //- Object name for name server registration.
        Foam::fileName objectName_;

        //- List of licenced host names.
        Foam::HashPtrTable<FoamXServer::HostDescriptor, Foam::string> hosts_;

        //- Process controls
        FoamX::ProcessControl procControl_;


    // Protected member functions.

        void getServerReference
        (
            const char *hostName,
            const Foam::fileName& objectKey,
            FoamXServer::CaseBrowser::ICaseBrowser_out caseBrowser
        );


public:

    // Constructors

        IHostBrowserImpl(Orb& orb, const Foam::stringList& args);


    // Destructor

        virtual ~IHostBrowserImpl();


    // Member Functions

        // Implemented IHostBrowser methods.

            virtual FoamXServer::HostDescriptorList* hosts();
            virtual void refreshHostList();
            virtual CORBA::Boolean isHostAlive(const char* hostname);
            virtual void hostIsAlive(const char* hostname);
            virtual void hostIsDead(const char* hostname);
            virtual void validate();
            virtual void save();
            virtual void close();

            virtual void openCaseBrowser
            (
                const char* hostname
            );

            virtual CORBA::Boolean getCaseBrowserReference
            (
                const char* hostname,
                FoamXServer::CaseBrowser::ICaseBrowser_out browserObj
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FoamX

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
