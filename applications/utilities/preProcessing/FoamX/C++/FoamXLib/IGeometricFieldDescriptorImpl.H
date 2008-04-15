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
    Foam::IGeometricFieldDescriptorImpl

Description

SourceFiles
    IGeometricFieldDescriptorImpl.C

\*---------------------------------------------------------------------------*/

#ifndef IGeometricFieldDescriptorImpl_H
#define IGeometricFieldDescriptorImpl_H

// Foam header files.
#include "IDLList.H"
#include "dictionary.H"

// FoamX header files.
#include "FoamX.H"
#include "FoamXErrors.H"
#include "DictionaryWriter.H"
#include "ObjRefHashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace FoamX
{

// Forward declaration of classes.
class ITypeDescriptorImpl;
class IGeometryDescriptorImpl;

/*---------------------------------------------------------------------------*\
                  Class IGeometricFieldDescriptorImpl Declaration
\*---------------------------------------------------------------------------*/

class IGeometricFieldDescriptorImpl
:
    public POA_FoamXServer::CaseServer::IGeometricFieldDescriptor,
    public PortableServer::RefCountServantBase,
    public IDLList<IGeometricFieldDescriptorImpl>::link
{
    // Private data

        //- Field name (eg, "p").
        Foam::word fieldName_;

        //- Field description (eg, "Pressure").
        Foam::string fieldDescription_;

        //- Field dimensions.
        FoamXServer::DimensionSet dimensionSet_;

        //- Field type name (eg, "scalar").
        Foam::word fieldTypeName_;

        //- Reference to type descriptor.
        ITypeDescriptorImpl* typeDescriptor_;

        //- Reference to type descriptor for the field.
        ITypeDescriptorImpl* fieldTypeDescriptor_;

        //- Geometry type name (eg, "vol").
        Foam::word geometryTypeName_;

        //- Reference tp geometry descriptor for this field.
        IGeometryDescriptorImpl* geometryDescriptor_;


public:

    // Constructors

        IGeometricFieldDescriptorImpl(const char* fieldName);

        IGeometricFieldDescriptorImpl
        (
            const Foam::word& fieldName,
            const Foam::dictionary& fieldDict,
            const ObjRefHashTable<ITypeDescriptorImpl*>& foamTypes,
            const ObjRefHashTable<IGeometryDescriptorImpl*>& geometryTypes
        );


    // Destructor

        virtual ~IGeometricFieldDescriptorImpl();


    // Member Functions
 
        // Foam access functions

            //- Keyword used by Dictionary
            const word& keyword() const
            {
                return fieldName_;
            }


        // FoamServer access functions

            virtual FoamXServer::ITypeDescriptor_ptr typeDescriptor();
            virtual FoamXServer::ITypeDescriptor_ptr fieldTypeDescriptor();

            virtual FoamXServer::CaseServer::IGeometryDescriptor_ptr
                geometryDescriptor();

            virtual char* name();
            virtual void name(const char*);

            virtual char* description();
            virtual void description(const char*);

            virtual char* fieldTypeName();
            virtual void fieldTypeName(const char*);

            virtual char* geometryTypeName();
            virtual void geometryTypeName(const char*);

            virtual FoamXServer::DimensionSet dimensions();
            virtual void dimensions
            (
                const FoamXServer::DimensionSet&
            );


        //- Persist to file.
        void save(DictionaryWriter& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FoamX

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
