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

#include "long.H"
#include "string.H"
#include "fileName.H"
#include "dictionary.H"
#include "LogEntry.H"

// Project header files.
#include "FoamX.H"
#include "FoamXErrors.H"

// Namespaces
#include "FoamXNameSpaces.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(FoamX::FoamXError, 0);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::FoamXError::FoamXError
(
    FoamXServer::ErrorCode errCode,
    const string& message,
    const char* function,
    const char* fName, 
    label lineNo
)
:
    FoamXServer::FoamXError
    (
        errCode,
        message.c_str(),
        function,
        fName, 
        lineNo
    )
{
    static const char* functionName =
        "FoamX::FoamXError::FoamXError"
        "(FoamXServer::ErrorCode errCode, const string& message,"
        "const char* functionName, const char* fName, label lineNo)";

    LogEntry log(functionName, __FILE__, __LINE__);
    log << *this << endl;
}

FoamX::FoamXError::FoamXError(const FoamXServer::FoamXError& fErr)
:
    FoamXServer::FoamXError(fErr)
{
    static const char* functionName =
        "FoamX::FoamXError::FoamXError(const FoamXServer::FoamXError& fErr)";

    LogEntry log(functionName, __FILE__, __LINE__);
    log << *this << endl;
}

FoamX::FoamXError::FoamXError(const error& fErr)
:
    FoamXServer::FoamXError
    (
        FoamXServer::E_FOAM,
        fErr.message().c_str(),
        fErr.functionName().c_str(),
        fErr.sourceFileName().c_str(),
        fErr.sourceFileLineNumber()
    )
{
    static const char* functionName =
        "FoamX::FoamXError::FoamXError(const error& fErr)";

    LogEntry log(functionName, __FILE__, __LINE__);
    log << *this << endl;
}


Foam::Ostream& FoamX::operator<<
(
    Ostream& os,
    const FoamXServer::FoamXError& fxErr
)
{
    os  << "FoamXError " << string(fxErr.errorMessage) << nl
        << "In function " << string(fxErr.methodName) << nl
        << "in file " << fileName(fxErr.fileName)
        << " at line " << fxErr.lineNo;

    return os;
}


FoamX::FoamXSYSError::FoamXSYSError
(
    FoamXServer::ErrorCode errCode,
    const string& message,
    const string& hostName,
    const char* function,
    const char* fName, 
    label lineNo
)
:
    FoamXServer::FoamXSYSError
    (
        errCode,
        message.c_str(),
        hostName.c_str(),
        function,
        fName, 
        lineNo
    )
{
    static const char* functionName =
        "FoamX::FoamXError::FoamXSYSError"
        "(FoamXServer::ErrorCode errCode, const string& messsage, "
        "const string& hostName, "
        "const char* functionName, const char* fName, label lineNo)";

    LogEntry log(functionName, __FILE__, __LINE__);
    log << *this << endl;
}

FoamX::FoamXSYSError::FoamXSYSError(const FoamXServer::FoamXSYSError& sysErr)
:
    FoamXServer::FoamXSYSError(sysErr)
{
    static const char* functionName =
        "FoamX::FoamXError::FoamXError"
        "(const FoamXServer::FoamXSYSError& sysErr)";

    LogEntry log(functionName, __FILE__, __LINE__);
    log << *this << endl;
}

Foam::Ostream& FoamX::operator<<
(
    Ostream& os,
    const FoamXServer::FoamXSYSError& sysErr
)
{
    os  << "FoamXSYSError " << string(sysErr.errorMessage) << nl
        << "Problem with machine " << string(sysErr.hostName) << nl
        << "In function " << string(sysErr.methodName) << nl
        << "in file " << fileName(sysErr.fileName)
        << " at line " << sysErr.lineNo;

    return os;
}


FoamX::FoamXIOError::FoamXIOError
(
    const string& message,
    const string& ioFileName,
    label ioStartLineNumber,
    label ioEndLineNumber,
    const char* function,
    const char* fName, 
    label lineNo
)
:
    FoamXServer::FoamXIOError
    (
        message.c_str(),
        ioFileName.c_str(),
        ioStartLineNumber,
        ioEndLineNumber,
        function,
        fName, 
        lineNo
    )
{
    static const char* functionName =
        "FoamX::FoamXError::FoamXError"
        "(const string& message, const string& ioFileName, "
        "label ioStartLineNumber, label ioEndLineNumber, "
        "const char* functionName, const char* fName, label lineNo)";

    LogEntry log(functionName, __FILE__, __LINE__);
    log << *this << endl;
}

FoamX::FoamXIOError::FoamXIOError(const FoamXServer::FoamXIOError& fErr)
:
    FoamXServer::FoamXIOError(fErr)
{
    static const char* functionName =
        "FoamX::FoamXIOError::FoamXIOError"
        "(const FoamXServer::FoamXIOError& fErr)";

    LogEntry log(functionName, __FILE__, __LINE__);
    log << *this << endl;
}

FoamX::FoamXIOError::FoamXIOError(const IOerror& fIOErr)
:
    FoamXServer::FoamXIOError
    (
        fIOErr.message().c_str(),
        fIOErr.ioFileName().c_str(),
        fIOErr.ioStartLineNumber(),
        fIOErr.ioEndLineNumber(),
        fIOErr.functionName().c_str(),
        fIOErr.sourceFileName().c_str(),
        fIOErr.sourceFileLineNumber()
    )
{
    static const char* functionName =
        "FoamX::foamXError(const IOerror& fIOErr)";

    LogEntry log(functionName, __FILE__, __LINE__);
    log << *this << endl;
}


Foam::Ostream& FoamX::operator<<
(
    Ostream& os,
    const FoamXServer::FoamXIOError& fIOErr
)
{
    os  << "FoamXIOError " << string(fIOErr.errorMessage) << nl
        << "File " << fileName(fIOErr.ioFileName)
        << " starting at line " << fIOErr.ioStartLineNo
        << " ending at line " << fIOErr.ioEndLineNo << nl
        << "In function " << string(fIOErr.methodName) << nl
        << "in file " << fileName(fIOErr.fileName)
        << " at line " << fIOErr.lineNo;

    return os;
}



Foam::dictionary FoamX::dict(const error& fErr)
{
    return fErr;
}

Foam::dictionary FoamX::dict(const IOerror& fIOErr)
{
    return fIOErr;
}

Foam::dictionary FoamX::dict(const FoamXServer::FoamXError& fxErr)
{
    dictionary fxErrDict;
    fxErrDict.add("type", word("FoamXServer::FoamXError"));
    fxErrDict.add("errorCode", label(fxErr.errorCode));
    fxErrDict.add("message", string(fxErr.errorMessage));
    fxErrDict.add("function", string(fxErr.methodName));
    fxErrDict.add("sourceFile", fileName(fxErr.fileName));
    fxErrDict.add("sourceFileLineNumber", fxErr.lineNo);

    return fxErrDict;
}

Foam::dictionary FoamX::dict(const FoamXServer::FoamXIOError& fIOErr)
{
    dictionary fIOErrDict;
    fIOErrDict.add("type", word("FoamXServer::FoamXIOError"));
    fIOErrDict.add("message", string(fIOErr.errorMessage));
    fIOErrDict.add("ioFileName", fileName(fIOErr.ioFileName));
    fIOErrDict.add("ioStartLineNumber", fIOErr.ioStartLineNo);
    fIOErrDict.add("ioEndLineNumber", fIOErr.ioEndLineNo);
    fIOErrDict.add("function", string(fIOErr.methodName));
    fIOErrDict.add("sourceFile", fileName(fIOErr.fileName));
    fIOErrDict.add("sourceFileLineNumber", fIOErr.lineNo);

    return fIOErrDict;
}

Foam::dictionary FoamX::dict(const CORBA::COMM_FAILURE& ex)
{
    dictionary exDict;
    exDict.add("type", word("CORBA::COMM_FAILURE"));
    return exDict;
}

Foam::dictionary FoamX::dict(const CORBA::SystemException& ex)
{
    dictionary exDict;
    exDict.add("type", word("CORBA::SystemException"));
    return exDict;
}

Foam::dictionary FoamX::dict(const CORBA::Exception& ex)
{
    dictionary exDict;
    exDict.add("type", word("CORBA::Exception"));
    return exDict;
}

Foam::dictionary FoamX::dict(const FoamX::systemError& ex)
{
    dictionary exDict;
    exDict.add("type", word("::systemError"));
    return exDict;
}


void FoamX::reThrow(const dictionary& errorDict)
{
    word errorType(errorDict.lookup("type"));

    if (errorType == "error")
    {
        throw FoamXError(error(errorDict));
    }
    else if (errorType == "IOerror")
    {
        throw FoamXIOError(IOerror(errorDict));
    }
    else if (errorType == "FoamXServer::FoamXError")
    {
        throw FoamXError
        (
            FoamXServer::ErrorCode
            (
                readLabel(errorDict.lookup("errorCode"))
            ),
            string(errorDict.lookup("message")),
            string(errorDict.lookup("function")).c_str(),
            string(errorDict.lookup("sourceFile")).c_str(),
            readLabel(errorDict.lookup("sourceFileLineNumber"))
        );
    }
    else if (errorType == "FoamXServer::FoamXIOError")
    {
        throw FoamXIOError
        (
            string(errorDict.lookup("message")),
            string(errorDict.lookup("ioFileName")),
            readLabel(errorDict.lookup("ioStartLineNumber")),
            readLabel(errorDict.lookup("ioEndLineNumber")),
            string(errorDict.lookup("function")).c_str(),
            string(errorDict.lookup("sourceFile")).c_str(),
            readLabel(errorDict.lookup("sourceFileLineNumber"))
        );
    }
    else
    {
        throw FoamXError
        (
            FoamXServer::E_UNEXPECTED,
            "Unexpected error: " + errorType,
            "FoamX::reThrow(const dictionary& errorDict)",
            __FILE__, __LINE__
        );
    }
}


// ************************************************************************* //
