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
    FoamXError

Description

SourceFiles
    FoamXErrors.C

\*---------------------------------------------------------------------------*/

#ifndef FoamXError_H
#define FoamXError_H

// Foam header files.
#include "error.H"
#include "dictionaryEntry.H"
#include "OSstream.H"
#include "OFstream.H"
#include "className.H"

// Project header files.
#include "FoamX.H"
#include "FoamXErrors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace FoamX
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class systemError
{};


class FoamXError
:
    public FoamXServer::FoamXError
{
public:

    // Declare name of the class and it's debug switch
    ClassName("FoamXError");

    FoamXError
    (
        FoamXServer::ErrorCode,
        const Foam::string& message,
        const char* functionName,
        const char* sourceName,
        Foam::label sourceLineNo
    );

    FoamXError(const FoamXServer::FoamXError&);

    FoamXError(const Foam::error&);
};

Foam::Ostream& operator<<
(
    Foam::Ostream& os,
    const FoamXServer::FoamXError& fxErr
);


class FoamXSYSError
:
    public FoamXServer::FoamXSYSError
{
public:

    FoamXSYSError
    (
        FoamXServer::ErrorCode,
        const Foam::string& message,
        const Foam::string& hostName,
        const char* functionName,
        const char* sourceName,
        Foam::label sourceLineNo
    );

    FoamXSYSError(const FoamXServer::FoamXSYSError&);

    FoamXSYSError(const Foam::error&);
};

Foam::Ostream& operator<<
(
    Foam::Ostream& os,
    const FoamXServer::FoamXSYSError& fIOErr
);


class FoamXIOError
:
    public FoamXServer::FoamXIOError
{
public:

    FoamXIOError
    (
        const Foam::string& message,
        const Foam::string& ioFileName,
        Foam::label ioStartLineNumber,
        Foam::label ioEndLineNumber,
        const char* functionName,
        const char* sourceName,
        Foam::label sourceLineNo
    );

    FoamXIOError(const FoamXServer::FoamXIOError&);

    FoamXIOError(const Foam::IOerror&);
};

Foam::Ostream& operator<<
(
    Foam::Ostream& os,
    const FoamXServer::FoamXIOError& fIOErr
);


Foam::dictionary dict(const Foam::error& fErr);
Foam::dictionary dict(const Foam::IOerror& fIOErr);
Foam::dictionary dict(const FoamXServer::FoamXError& fxErr);
Foam::dictionary dict(const FoamXServer::FoamXIOError& fxIOErr);
Foam::dictionary dict(const CORBA::COMM_FAILURE& ex);
Foam::dictionary dict(const CORBA::SystemException& ex);
Foam::dictionary dict(const CORBA::Exception& ex);
Foam::dictionary dict(const FoamX::systemError& ex);


void reThrow(const Foam::dictionary& errorDict);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template <class Exception>
inline void writeException(const Exception& ex, const Foam::fileName& fName)
{
    if (fName.size())
    {
        Foam::OFstream ofStream(fName);

        if (ofStream.good())
        {
            ofStream
                << Foam::dictionaryEntry
                   (
                       "FatalError",
                       Foam::dictionary::null,
                       dict(ex)
                   )
                << Foam::endl;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End of namespace FoamX

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define CATCH_ALL(functionName)                                              \
                                                                             \
    catch (Foam::IOerror& fIOErr)                                            \
    {                                                                        \
        /* Convert the Foam::IOerror into a FoamX::IOerror error             \
           and re-throw */                                                   \
        throw FoamX::FoamXIOError(fIOErr);                                   \
    }                                                                        \
    catch (Foam::error& fErr)                                                \
    {                                                                        \
        /* Convert the Foam error into a FoamX error and re-throw */         \
        throw FoamX::FoamXError(fErr);                                       \
    }                                                                        \
    catch (FoamX::FoamXError& fxErr)                                         \
    {                                                                        \
        /* Bounce exception up to client */                                  \
        throw fxErr;                                                         \
    }                                                                        \
    catch (FoamXServer::FoamXError& fxErr)                                   \
    {                                                                        \
        /* Convert into FoamX error and bounce up to client */               \
        throw FoamX::FoamXError(fxErr);                                      \
    }                                                                        \
    catch (FoamX::FoamXSYSError& sysErr)                                     \
    {                                                                        \
        /* Bounce exception up to client */                                  \
        throw sysErr;                                                        \
    }                                                                        \
    catch (FoamXServer::FoamXSYSError& sysErr)                               \
    {                                                                        \
        /*  Convert into FoamX error and bounce up to client */              \
        throw FoamX::FoamXSYSError(sysErr);                                  \
    }                                                                        \
    catch (FoamX::FoamXIOError& fIOErr)                                      \
    {                                                                        \
        /* Bounce exception up to client */                                  \
        throw fIOErr;                                                        \
    }                                                                        \
    catch (FoamXServer::FoamXIOError& fIOErr)                                \
    {                                                                        \
        /* Convert into FoamX error and bounce up to client */               \
        throw FoamX::FoamXIOError(fIOErr);                                   \
    }                                                                        \
    catch (...)                                                              \
    {                                                                        \
        throw FoamX::FoamXError                                              \
        (                                                                    \
            FoamXServer::E_UNEXPECTED,                                       \
            "Unexpected error.",                                             \
            functionName,                                                    \
            __FILE__, __LINE__                                               \
        );                                                                   \
    }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
