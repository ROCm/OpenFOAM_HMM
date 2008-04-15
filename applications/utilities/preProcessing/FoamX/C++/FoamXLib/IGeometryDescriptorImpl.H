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
    Foam::IPatchDescriptorImpl

Description

SourceFiles
    IGeometryDescriptorImpl.C

\*---------------------------------------------------------------------------*/

#ifndef IGeometryDescriptorImpl_H
#define IGeometryDescriptorImpl_H

// Foam header files.
#include "word.H"
#include "string.H"
#include "dictionary.H"

// FoamX header files.
#include "FoamX.H"
#include "FoamXErrors.H"
#include "DictionaryWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace FoamX
{

/*---------------------------------------------------------------------------*\
                 Class IGeometryDescriptorImpl Declaration
\*---------------------------------------------------------------------------*/

class IGeometryDescriptorImpl
:
    public POA_FoamXServer::CaseServer::IGeometryDescriptor,
    public PortableServer::RefCountServantBase
{
    // Private data

        Foam::word name_;
        Foam::string displayName_;
        Foam::string description_;


public:

    // Constructors

        IGeometryDescriptorImpl(const word& geometryName);

        IGeometryDescriptorImpl
        (
            const word& geometryName,
            const dictionary& outerDict
        );


    // Destructor

        virtual ~IGeometryDescriptorImpl();


    // Member Functions
  
        // Properties common to all types.

            virtual char* name();
            virtual void name(const char*);
            virtual char* displayName();
            virtual void displayName(const char*);
            virtual char* description();
            virtual void description(const char*);


        //- Persist to file.
        void save(DictionaryWriter& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FoamX

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
