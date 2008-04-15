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
    NameServer.C

\*---------------------------------------------------------------------------*/

#ifndef NameServer_H
#define NameServer_H

// System header files.
#include <coss/CosNaming.h>

// Foam header files.
#include "dictionary.H"

// FoamX header files.
#include "FoamX.H"
#include "FoamXErrors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace FoamX
{

/*---------------------------------------------------------------------------*\
                           Class NameServer Declaration
\*---------------------------------------------------------------------------*/

class NameServer
{
    // Private data

        bool connected_;
        CosNaming::NamingContext_var rootContext_;


    // Private Member Functions

        void createNameFromString(const std::string&, CosNaming::Name& name);


public:

    // Constructors

        //- Construct null
        NameServer();

        //- Construct connecting to given ORB
        NameServer(CORBA::ORB_ptr pOrb);


    // Destructor

        virtual ~NameServer();


    // Member Functions

        // Connection methods.

            void connect(CORBA::ORB_ptr pOrb);

            void disconnect();

            bool isConnected() const
            {
                return connected_;
            }


        // Context methods.

            // Create single context
            void createContext
            (
                const Foam::fileName& path,
                bool failAlreadyBound = false
            );

            // Create all contexts (a la 'mkdir -p')
            void createContexts
            (
                const Foam::fileName& path,
                bool failAlreadyBound = false
            );

            void removeContext(const Foam::fileName& path);

            // Remove all contexts up to (but not including) root
            void removeContexts
            (
                const Foam::fileName& root,
                const Foam::fileName& path
            );


            bool isContextBound(const Foam::fileName& path);


        // Object binding methods.

            void bindObject
            (
                const Foam::fileName& name,
                CORBA::Object_ptr pObject,
                bool rebind = false
            );

            void unbindObject(const Foam::fileName& name);

            bool isObjectBound(const Foam::fileName& name);


        // Object reference resolution methods.

            CORBA::Object_ptr resolve(const Foam::fileName& name);


        //- Template function to return object reference of the template type
        //  from the given context name.
        //  For Example :-
        //  FoamXServer::HostBrowser::IHostBrowser_var ref =
        //      fxNameServer.resolve
        //      <FoamXServer::HostBrowser::IHostBrowser>
        //      ("HostBrowser.FoamX");
        template<class T>
        typename T::_ptr_type resolve(const Foam::fileName& name)
        {
            CORBA::Object_var objRef = resolve(name);

            // Narrow to the specified type.
            typename T::_ptr_type objTypeRef = T::_narrow(objRef);
            if (CORBA::is_nil(objTypeRef))
            {
                throw FoamXError
                (
                    FoamXServer::E_FAIL,
                    "Invalid object type.",
                    "NameServer::resolve(const Foam::fileName& name)",
                    __FILE__, __LINE__
                );
            }

            return  T::_duplicate(objTypeRef);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FoamX

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
