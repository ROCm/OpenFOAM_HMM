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
#include "emptyPolyPatch.H"
#include "polyPatch.H"

// FoamX header files.
#include "FoamX.H"
#include "LogEntry.H"
#include "PatchProperties.H"
#include "Paths.H"

// Namespaces
#include "FoamXNameSpaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::PatchProperties::PatchProperties(const word& patchName)
:
    patchName_(patchName),
    patchType_(emptyPolyPatch::typeName),
    physicalType_(emptyPolyPatch::typeName),
    modified_(false)
{
    static const char* functionName =
        "FoamX::PatchProperties::PatchProperties(const char* patchName)";

    LogEntry log(functionName, __FILE__, __LINE__);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

FoamX::PatchProperties::~PatchProperties()
{
    static const char* functionName =
        "FoamX::PatchProperties::~PatchProperties()";

    LogEntry log(functionName, __FILE__, __LINE__);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::word& FoamX::PatchProperties::patchName() const
{
    static const char* functionName =
        "FoamX::PatchProperties::patchName()";

    LogEntry log(functionName, __FILE__, __LINE__);

    return patchName_;
}

void FoamX::PatchProperties::patchName(const word& newName)
{
    static const char* functionName =
        "FoamX::PatchProperties::patchName(const char* newName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    patchName_ = newName;

    modified_ = true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::word& FoamX::PatchProperties::patchType() const
{
    static const char* functionName =
        "FoamX::PatchProperties::patchType()";

    LogEntry log(functionName, __FILE__, __LINE__);

    return patchType_;
}

void FoamX::PatchProperties::patchType(const word& newPatchType)
{
    static const char* functionName =
        "FoamX::PatchProperties::type(const char* newPatchType)";

    LogEntry log(functionName, __FILE__, __LINE__);

    patchType_ = newPatchType;

    modified_ = true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::word& FoamX::PatchProperties::physicalType() const
{
    static const char* functionName =
        "FoamX::PatchProperties::physicalType()";

    LogEntry log(functionName, __FILE__, __LINE__);

    return physicalType_;
}

void FoamX::PatchProperties::physicalType(const word& newName)
{
    static const char* functionName =
        "FoamX::PatchProperties::physicalType(const char* newName)";

    LogEntry log(functionName, __FILE__, __LINE__);

    physicalType_ = newName;

    modified_ = true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::label FoamX::PatchProperties::startFace() const
{
    static const char* functionName =
        "FoamX::PatchProperties::startFace()";

    LogEntry log(functionName, __FILE__, __LINE__);

    return startFace_;
}

void FoamX::PatchProperties::startFace(label n)
{
    static const char* functionName =
        "FoamX::PatchProperties::startFace(label n)";

    LogEntry log(functionName, __FILE__, __LINE__);

    startFace_ = n;

    modified_ = true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::label FoamX::PatchProperties::nFaces() const
{
    static const char* functionName =
        "FoamX::PatchProperties::nFaces()";

    LogEntry log(functionName, __FILE__, __LINE__);

    return nFaces_;
}

void FoamX::PatchProperties::nFaces(label n)
{
    static const char* functionName =
        "FoamX::PatchProperties::nFaces(label n)";

    LogEntry log(functionName, __FILE__, __LINE__);

    nFaces_ = n;

    modified_ = true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool FoamX::PatchProperties::modified() const
{
    static const char* functionName =
        "FoamX::PatchProperties::modified()";

    LogEntry log(functionName, __FILE__, __LINE__);

    return modified_;
}

void FoamX::PatchProperties::modified(bool mod)
{
    static const char* functionName =
        "FoamX::PatchProperties::modified(bool)";

    LogEntry log(functionName, __FILE__, __LINE__);

    modified_ = mod;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::PatchProperties::load(const dictionary& patchDict)
{
    static const char* functionName =
        "FoamX::PatchProperties::load(const dictionary& patchDict)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        log << "Reading " << patchDict.name()
            << " start line " << patchDict.startLineNumber()
            << " end line " << patchDict.endLineNumber() << endl;

        // Check that the dictionary has the minimum required information.
        if
        (
            !patchDict.found("type")
         || !patchDict.found("startFace")
         || !patchDict.found("nFaces")
        )
        {
            throw FoamXError
            (
                E_FAIL,
                "Malformed patch dictionary '" + patchName_ + "'.",
                functionName,
                __FILE__, __LINE__
            );
        }

        // Get patch data.
        patchDict.lookup("type") >> patchType_;

        if (patchDict.found("geometricType"))
        {
            patchDict.lookup("geometricType") >> patchType_;
        }

        patchDict.lookup("startFace") >> startFace_;
        patchDict.lookup("nFaces") >> nFaces_;

        // If the type is a constraint type set the physicalType
        // to be the same else lookup the physical type if it is there.
        if (polyPatch::constraintType(patchType_))
        {
            physicalType_ = patchType_;
        }
        else if (patchDict.found("physicalType"))
        {
            patchDict.lookup("physicalType") >> physicalType_;
        }

        modified_ = false;
    }
    catch (FoamXError& ex)
    {
        // Bounce exception up to client.
        throw ex;
    }
    catch (...)
    {
        throw FoamXError
        (
            E_UNEXPECTED,
            "Unexpected error.",
            functionName,
            __FILE__, __LINE__
        );
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void FoamX::PatchProperties::save(FoamX::DictionaryWriter& dictWriter)
{
    static const char* functionName =
        "FoamX::PatchProperties::save(FoamX::DictionaryWriter& dictWriter)";

    LogEntry log(functionName, __FILE__, __LINE__);

    try
    {
        // Save the patch information to a sub-dictionary.
        dictWriter.startSubDict(patchName_);

        dictWriter.writeEntry("type", patchType_);
        dictWriter.writeEntry("physicalType", physicalType_);
        dictWriter.writeEntry("startFace", startFace_  );
        dictWriter.writeEntry("nFaces", nFaces_);

        dictWriter.endSubDict();

        modified_ = false;
    }
    catch (FoamXError& ex)
    {
        // Bounce exception up to client.
        throw ex;
    }
    catch (...)
    {
        throw FoamXError
        (
            E_UNEXPECTED,
            "Unexpected error.",
            functionName,
            __FILE__, __LINE__
        );
    }
}


// ************************************************************************* //
