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
#include "FoamXString.H"

// Namespaces
#include "FoamXNameSpaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::FoamXString::FoamXString()
{
    tokenType_ = token::STRING;
    foamXType_ = FoamXServer::Type_String;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::FoamXString::FoamXString(const FoamX::FoamXString& str)
:
    String_var(str)
{
    tokenType_ = str.tokenType_;
    foamXType_ = str.foamXType_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::FoamXString::FoamXString(const Foam::word& str)
:
    String_var(str.c_str())
{
    tokenType_ = token::WORD;
    foamXType_ = FoamXServer::Type_Word;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::FoamXString::FoamXString(const Foam::string& str)
:
    String_var(str.c_str())
{
    tokenType_ = token::STRING;
    foamXType_ = FoamXServer::Type_String;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::FoamXString::FoamXString(const char* str)
:
    String_var(str)
{
    tokenType_ = token::STRING;
    foamXType_ = FoamXServer::Type_String;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::FoamXString::FoamXString(Foam::Istream& is)
{
    read(is);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::FoamXString::~FoamXString()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::FoamXString& FoamX::FoamXString::operator=(const FoamXString& str)
{
    tokenType_ = str.tokenType_;
    foamXType_ = str.foamXType_;

    // Copy string value. Need to call base classes assignment operator.
    CORBA::String_var::operator=(str);
    return *this;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::FoamXString& FoamX::FoamXString::operator=(const Foam::string& str)
{
    tokenType_ = token::STRING;
    foamXType_ = FoamXServer::Type_String;

    // Copy string value. Need to call base classes assignment operator.
    CORBA::String_var::operator=(str.c_str());
    return *this;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::FoamXString& FoamX::FoamXString::operator=(const Foam::word& str)
{
    tokenType_ = token::WORD;
    foamXType_ = FoamXServer::Type_Word;

    // Copy string value. Need to call base classes assignment operator.
    CORBA::String_var::operator=(str.c_str());
    return *this;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::FoamXString& FoamX::FoamXString::operator=(const char* str)
{
    // Do not reset the type.

    // Copy string value. Need to call base classes assignment operator.
    CORBA::String_var::operator=(str);
    return *this;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::FoamXString::read(Foam::Istream& is)
{
    // Construct token from Istream.
    token tok(is);

    // Store the token type.
    tokenType_ = tok.type();

    if (tokenType_ == token::WORD)
    {
        CORBA::String_var::operator=(tok.wordToken().c_str());
        foamXType_ = FoamXServer::Type_Word;
    }
    else if (tokenType_ == token::STRING)
    {
        CORBA::String_var::operator=(tok.stringToken().c_str());
        foamXType_ = FoamXServer::Type_String;
    }
    else
    {
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::FoamXString::write(Foam::Ostream& os) const
{
    // Convert this string to the appropriate token type and output.

    if (tokenType_ == token::WORD)
    {
        os << Foam::word(this->in());       // In method returns const char*.
    }
    else if (tokenType_ == token::STRING)
    {
        os << Foam::string(this->in());     // In method returns const char*.
    }
    else
    {
        WarningIn("FoamX::FoamXString::write(Foam::Ostream& os) const")
            << "Expected word or string token."
            << Foam::endl;
    }
}


// ************************************************************************* //
