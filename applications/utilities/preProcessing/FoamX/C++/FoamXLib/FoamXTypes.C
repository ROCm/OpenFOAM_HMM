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

#include "FoamXTypes.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::wordList FoamX::FoamXTypes::typeNames_;
Foam::HashTable<FoamXServer::FoamXType> FoamX::FoamXTypes::types_;
Foam::word FoamX::FoamXTypes::invalidTypeName_("invalidType");

//- Construct dummy FoamXTypes to initialise static data
FoamX::FoamXTypes FoamXTypes_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct null
FoamX::FoamXTypes::FoamXTypes()
{
    typeNames_.setSize(FoamXServer::Type_Field + 1);

    typeNames_[FoamXServer::Type_Boolean] = "boolean";
    typeNames_[FoamXServer::Type_Label] = "label";
    typeNames_[FoamXServer::Type_Scalar] = "scalar";
    typeNames_[FoamXServer::Type_Char] = "char";
    typeNames_[FoamXServer::Type_Word] = "word";
    typeNames_[FoamXServer::Type_String] = "string";
    typeNames_[FoamXServer::Type_RootDir] = "rootDir";
    typeNames_[FoamXServer::Type_RootAndCase] = "rootAndCase";
    typeNames_[FoamXServer::Type_CaseName] = "caseName";
    typeNames_[FoamXServer::Type_HostName] = "hostName";
    typeNames_[FoamXServer::Type_File] = "file";
    typeNames_[FoamXServer::Type_Directory] = "directory";
    typeNames_[FoamXServer::Type_Time] = "time";
    typeNames_[FoamXServer::Type_DimensionSet] = "dimensionSet";
    typeNames_[FoamXServer::Type_FixedList] = "fixedList";
    typeNames_[FoamXServer::Type_List] = "list";
    typeNames_[FoamXServer::Type_Dictionary] = "dictionary";
    typeNames_[FoamXServer::Type_Selection] = "selection";
    typeNames_[FoamXServer::Type_Compound] = "compound";
    typeNames_[FoamXServer::Type_Field] = "field";

    forAll(typeNames_, i)
    {
        types_.insert(typeNames_[i], FoamXServer::FoamXType(i));
    }
}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

bool FoamX::FoamXTypes::found(const Foam::word& typeName)
{
    return types_.found(typeName);
}

const Foam::word& FoamX::FoamXTypes::typeName(FoamXServer::FoamXType fxType)
{
    if
    (
        fxType >= FoamXServer::Type_Boolean
     && fxType <= FoamXServer::Type_Field
    )
    {
        return typeNames_[fxType];
    }
    else
    {
        return invalidTypeName_;
    }
}

FoamXServer::FoamXType FoamX::FoamXTypes::lookupType(const Foam::word& typeName)
{
    return types_[typeName];
}

bool FoamX::FoamXTypes::isNumber(FoamXServer::FoamXType type)
{
    return
    (
        type == FoamXServer::Type_Label || type == FoamXServer::Type_Scalar
    );
}

bool FoamX::FoamXTypes::isPrimitive(FoamXServer::FoamXType type)
{
    return (type <= FoamXServer::Type_DimensionSet);
}

bool FoamX::FoamXTypes::isPrimitive(const Foam::word& typeName)
{
    if (types_.found(typeName))
    {
        return isPrimitive(*(types_.find(typeName)));
    }
    else
    {
        return false;
    }
}

bool FoamX::FoamXTypes::isCompound(FoamXServer::FoamXType type)
{
    return 
    (
        type >= FoamXServer::Type_FixedList
     && type <= FoamXServer::Type_Compound
    );
}

bool FoamX::FoamXTypes::isCompound(const Foam::word& typeName)
{
    if (types_.found(typeName))
    {
        return isCompound(*(types_.find(typeName)));
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
