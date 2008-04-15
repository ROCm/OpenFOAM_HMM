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

// Foam header files
#include "dimensionSet.H"

// FoamX header files.
#include "FoamXAny.H"
#include "FoamXTypes.H"
#include "DimensionSet.H"

// Namespaces
#include "FoamXNameSpaces.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
// Copy constructor.

FoamX::FoamXAny::FoamXAny(const FoamXAny& any)
:
    FoamXServer::FoamXAny()
{
    // Copy value and type via assignment operator.
    (*this) = any;
}

// Default type is Type_Undefined.
FoamX::FoamXAny::FoamXAny(FoamXServer::FoamXType t)
{
    type = t;
    setDefaultValue();
}

// Read constructor given type
FoamX::FoamXAny::FoamXAny(FoamXServer::FoamXType t, Foam::Istream& is)
{
    type = t;
    read(is);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

FoamX::FoamXAny::~FoamXAny()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void FoamX::FoamXAny::setType(FoamXServer::FoamXType t)
{
    type = t;
    setDefaultValue();
}


void FoamX::FoamXAny::setValue(const FoamXServer::FoamXAny& any)
{
    static const char* functionName =
        "FoamX::FoamXAny::SetValue(const FoamXServer::FoamXAny& any)";

    try
    {
        if (type != any.type)
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Type mismatch, expected " + FoamXTypes::typeName(type)
              + ", got " + FoamXTypes::typeName(any.type),
                functionName,
                __FILE__, __LINE__
            );
        }

        FoamXServer::FoamXAny::operator=(any);
    }
    CATCH_ALL(functionName);
}


void FoamX::FoamXAny::setMin()
{
    static const char* functionName =
        "FoamX::FoamXAny::setMin()";

    try
    {
        switch(type)
        {
            case Type_Label:
            {
                value <<= CORBA::Long(Foam::labelMin);
            }
            break;

            case Type_Scalar:
            {
                value <<= CORBA::Double(-VGREAT);
            }
            break;

            default:
                throw FoamXError
                (
                    E_UNEXPECTED,
                    "Invalid type for setMin: "
                    "expected a scalar or label, found "
                  + FoamXTypes::typeName(type),
                    functionName,
                    __FILE__, __LINE__
                );
        }
    }
    CATCH_ALL(functionName);
}


void FoamX::FoamXAny::setMax()
{
    static const char* functionName =
        "FoamX::FoamXAny::setMin()";

    try
    {
        switch(type)
        {
            case Type_Label:
            {
                value <<= CORBA::Long(Foam::labelMax);
            }
            break;

            case Type_Scalar:
            {
                value <<= CORBA::Double(VGREAT);
            }
            break;

            default:
                throw FoamXError
                (
                    E_UNEXPECTED,
                    "Invalid type for setMax: "
                    "expected a scalar or label, found "
                  + FoamXTypes::typeName(type),
                    functionName,
                    __FILE__, __LINE__
                );
        }
    }
    CATCH_ALL(functionName);
}


void FoamX::FoamXAny::read(Istream& is)
{
    static const char* functionName =
        "FoamX::FoamXAny::read(Istream& is)";

    try
    {
        switch(type)
        {
            case Type_Boolean:
            {
                token t(is);
                // Boolean -> tk_boolean any type.
                value <<= CORBA::Any::from_boolean(bool(t.labelToken()));
            }
            break;

            case Type_Char:
            case Type_Word:
            {
                token t(is);
                // Char, Word -> tk_string any type.
                value <<= t.wordToken().c_str();
            }
            break;

            case Type_String:
            case Type_RootDir:
            case Type_RootAndCase:
            case Type_CaseName:
            case Type_HostName:
            case Type_File:
            case Type_Directory:
            case Type_Time:
            {
                token t(is);
                // String -> tk_string any type.
                value <<= t.stringToken().c_str();
            }
            break;

            case Type_Label:
            {
                token t(is);
                // Label -> tk_long any type.
                value <<= CORBA::Long(t.labelToken());
            }
            break;

            case Type_Scalar:
            {
                token t(is);
                // Scalar -> tk_double any type.
                value <<= t.number();
            }
            break;

            case Type_DimensionSet:
            {
                // Construct a dimensionSet object and copy into
                // DimensionSet structure.
                DimensionSet fdm;
                fdm == dimensionSet(is);
                value <<= fdm;
            }
            break;

            case Type_FixedList:
            case Type_List:
            case Type_Dictionary:
            case Type_Compound:
            case Type_Selection:
            case Type_Undefined:
            case Type_Field:
            default:
                throw FoamXError
                (
                    E_INVALID_ARG,
                    "Invalid type " + FoamXTypes::typeName(type),
                    functionName,
                    __FILE__, __LINE__
                );
        }

        // Check stream status.
        is.check("FoamX::FoamXAny::read(Istream&)");
    }
    CATCH_ALL(functionName);
}


void FoamX::FoamXAny::write(Foam::Ostream& os) const
{
    static const char* functionName =
        "FoamX::FoamXAny::write(Foam::Ostream& os) const";

    try
    {
        if (!IsSet())
        {
            os << "<UnsetValue>";
            return;
        }

        // Write current value.
        switch(type)
        {
            case Type_Boolean:
            {
                // Extract into a boolean.
                CORBA::Boolean b = 0;
                CORBA::Any::to_boolean tb(b);
                if (!(value >>= tb))
                {
                    throw FoamXError
                    (
                        E_UNEXPECTED,
                        "Invalid type for Boolean output.",
                        functionName,
                        __FILE__, __LINE__
                    );
                }
                os << Foam::label(b);
            }
            break;

            case Type_Char:
            case Type_Word:
            {
                // Extract into a string.
                const char* str = NULL;
                if (!(value >>= str))
                {
                    throw FoamXError
                    (
                        E_UNEXPECTED,
                        "Invalid type for Char/Word output.",
                        functionName,
                        __FILE__, __LINE__
                    );
                }
                os << Foam::word(str);
            }
            break;

            case Type_String:
            case Type_RootDir:
            case Type_RootAndCase:
            case Type_CaseName:
            case Type_HostName:
            case Type_File:
            case Type_Time:
            {
                // Extract into a string.
                const char* str = NULL;
                if (!(value >>= str))
                {
                    throw FoamXError
                    (
                        E_UNEXPECTED,
                        "Invalid type for String output.",
                        functionName,
                        __FILE__, __LINE__
                    );
                }
                os << Foam::string(str);
            }
            break;

            case Type_Label:
            {
                // Extract into a long int.
                CORBA::Long lb = 0;
                if (!(value >>= lb))
                {
                    throw FoamXError
                    (
                        E_UNEXPECTED,
                        "Invalid type for Label output.",
                        functionName,
                        __FILE__, __LINE__
                    );
                }
                os << Foam::label(lb);
            }
            break;

            case Type_Scalar:
            {
                // Extract into a scalar.
                CORBA::Double sc = 0.0;
                if (!(value >>= sc))
                {
                    throw FoamXError
                    (
                        E_UNEXPECTED,
                        "Invalid type for Scalar output.",
                        functionName,
                        __FILE__, __LINE__
                    );
                }
                os << Foam::scalar(sc);
            }
            break;

            case Type_DimensionSet:
            {
                // Extract into dimension set.
                const DimensionSet* fdm = NULL;
                if (!(value >>= fdm))
                {
                    throw FoamXError
                    (
                        E_UNEXPECTED,
                        "Invalid type for DimensionSet output.",
                        functionName,
                        __FILE__, __LINE__
                    );
                }
                os  << token::BEGIN_SQR
                    <<        fdm->mass
                    << " " << fdm->length
                    << " " << fdm->time
                    << " " << fdm->temperature
                    << " " << fdm->moles
                    << " " << fdm->current
                    << " " << fdm->luminousIntensity
                    << token::END_SQR;
            }
            break;

            case Type_FixedList:
            case Type_List:
            case Type_Dictionary:
            case Type_Compound:
            case Type_Selection:
            case Type_Undefined:
            default:
                throw FoamXError
                (
                    E_UNEXPECTED,
                    "Invalid type.",
                    functionName,
                    __FILE__, __LINE__
                );
        }

        // Check stream status.
        os.check("FoamX::FoamXAny::write(Ostream&)");
    }
    CATCH_ALL(functionName);
}


void FoamX::FoamXAny::setDefaultValue()
{
    static const char* functionName =
        "FoamX::FoamXAny::setDefaultValue()";

    try
    {
        switch(type)
        {
            case Type_Boolean:
            {
                value <<= CORBA::Any::from_boolean(false);
            }
            break;

            case Type_Char:
            {
                value <<= static_cast<const char*>("?");
            }
            break;

            case Type_Word:
            {
                value <<= static_cast<const char*>("<>");
            }
            break;

            case Type_String:
            case Type_RootDir:
            case Type_RootAndCase:
            case Type_CaseName:
            case Type_HostName:
            case Type_File:
            case Type_Time:
            {
                value <<= static_cast<const char*>("");
            }
            break;

            case Type_Label:
            {
                value <<= CORBA::Long(0);
            }
            break;

            case Type_Scalar:
            {
                value <<= CORBA::Double(0.0);
            }
            break;

            case Type_DimensionSet:
            {
                DimensionSet fdm;
                fdm.mass              = 0;
                fdm.length            = 0;
                fdm.time              = 0;
                fdm.temperature       = 0;
                fdm.moles             = 0;
                fdm.current           = 0;
                fdm.luminousIntensity = 0;

                value <<= fdm;
            }
            break;

            default:
                // No default value required for any other types.
            break;
        }
    }
    CATCH_ALL(functionName);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void FoamX::FoamXAny::operator=(const FoamXServer::FoamXAny& any)
{
    static const char* functionName =
        "FoamX::FoamXAny::operator=(const FoamXServer::FoamXAny& any)";

    try
    {
        if (type != any.type)
        {
            throw FoamXError
            (
                E_INVALID_ARG,
                "Type mismatch, expected " + FoamXTypes::typeName(type)
              + ", got " + FoamXTypes::typeName(any.type),
                functionName,
                __FILE__, __LINE__
            );
        }

        FoamXServer::FoamXAny::operator=(any);
    }
    CATCH_ALL(functionName);
}


void FoamX::FoamXAny::operator=(const FoamXAny& any)
{
    // Check for assignment to self.
    if (this != & any)
    {
        // Copy type.
        type = any.type;

        // Copy value
        value = any.value;
    }
}


bool FoamX::FoamXAny::operator==(const FoamXServer::FoamXAny& any) const
{
    return (type == any.type && value == any.value);
}


bool FoamX::FoamXAny::operator!=(const FoamXServer::FoamXAny& any) const
{
    return !operator==(any);
}


void FoamX::operator>
(
    const FoamXServer::FoamXAny& any1,
    const FoamXServer::FoamXAny& any2
)
{
    static const char* functionName =
        "FoamX::operator>"
        "(const FoamXServer::FoamXAny& any1 const FoamXServer::FoamXAny& any2)";

    try
    {
        if (any1.type != any2.type)
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Incompatible types for operation >. "
                "Expected scalars or labels, found a "
              + FoamXTypes::typeName(any1.type) + " and a "
              + FoamXTypes::typeName(any2.type), 
                functionName,
                __FILE__, __LINE__
            );
        }

        switch(any1.type)
        {
            case Type_Label:
            {
                // Extract this into a long int.
                CORBA::Long lb1 = 0;
                if (!(any1.value >>= lb1))
                {
                    throw FoamXError
                    (
                        E_UNEXPECTED,
                        "Invalid type for Label output.",
                        functionName,
                        __FILE__, __LINE__
                    );
                }

                // Extract the any into a long int.
                CORBA::Long lb2 = 0;
                if (!(any2.value >>= lb2))
                {
                    throw FoamXError
                    (
                        E_UNEXPECTED,
                        "Invalid type for Label output.",
                        functionName,
                        __FILE__, __LINE__
                    );
                }

                if (lb1 > lb2)
                {
                    throw FoamXError
                    (
                        E_UNEXPECTED,
                        "Label value " + name(label(lb1))
                      + " is greater than the maximum " + name(label(lb2)),
                        functionName,
                        __FILE__, __LINE__
                    );
                }
            }
            break;

            case Type_Scalar:
            {
                // Extract this into a scalar.
                CORBA::Double sc1 = 0.0;
                if (!(any1.value >>= sc1))
                {
                    throw FoamXError
                    (
                        E_UNEXPECTED,
                        "Invalid type for Scalar output.",
                        functionName,
                        __FILE__, __LINE__
                    );
                }

                // Extract the any into a scalar.
                CORBA::Double sc2 = 0.0;
                if (!(any2.value >>= sc2))
                {
                    throw FoamXError
                    (
                        E_UNEXPECTED,
                        "Invalid type for Scalar output.",
                        functionName,
                        __FILE__, __LINE__
                    );
                }

                if (sc1 > sc2)
                {
                    throw FoamXError
                    (
                        E_UNEXPECTED,
                        "Label value " + name(scalar(sc1))
                      + " is greater than the maximum " + name(scalar(sc2)),
                        functionName,
                        __FILE__, __LINE__
                    );
                }
            }
            break;

            default:
                throw FoamXError
                (
                    E_UNEXPECTED,
                    "Invalid type for operation >.  "
                    "Expected scalars or labels, found "
                  + FoamXTypes::typeName(any1.type),
                    functionName,
                    __FILE__, __LINE__
                );
        }
    }
    CATCH_ALL(functionName);
}


void FoamX::operator<
(
    const FoamXServer::FoamXAny& any1,
    const FoamXServer::FoamXAny& any2
)
{
    static const char* functionName =
        "FoamX::operator<"
        "(const FoamXServer::FoamXAny& any1 const FoamXServer::FoamXAny& any2)";

    try
    {
        if (any1.type != any2.type)
        {
            throw FoamXError
            (
                E_UNEXPECTED,
                "Incompatible types for operation <. "
                "Expected scalars or labels, found a "
              + FoamXTypes::typeName(any1.type) + " and a "
              + FoamXTypes::typeName(any2.type), 
                functionName,
                __FILE__, __LINE__
            );
        }

        switch(any1.type)
        {
            case Type_Label:
            {
                // Extract this into a long int.
                CORBA::Long lb1 = 0;
                if (!(any1.value >>= lb1))
                {
                    throw FoamXError
                    (
                        E_UNEXPECTED,
                        "Invalid type for Label output.",
                        functionName,
                        __FILE__, __LINE__
                    );
                }

                // Extract the any into a long int.
                CORBA::Long lb2 = 0;
                if (!(any2.value >>= lb2))
                {
                    throw FoamXError
                    (
                        E_UNEXPECTED,
                        "Invalid type for Label output.",
                        functionName,
                        __FILE__, __LINE__
                    );
                }

                if (lb1 < lb2)
                {
                    throw FoamXError
                    (
                        E_UNEXPECTED,
                        "Label value " + name(label(lb1))
                      + " is greater than the maximum " + name(label(lb2)),
                        functionName,
                        __FILE__, __LINE__
                    );
                }
            }
            break;

            case Type_Scalar:
            {
                // Extract this into a scalar.
                CORBA::Double sc1 = 0.0;
                if (!(any1.value >>= sc1))
                {
                    throw FoamXError
                    (
                        E_UNEXPECTED,
                        "Invalid type for Scalar output.",
                        functionName,
                        __FILE__, __LINE__
                    );
                }

                // Extract the any into a scalar.
                CORBA::Double sc2 = 0.0;
                if (!(any2.value >>= sc2))
                {
                    throw FoamXError
                    (
                        E_UNEXPECTED,
                        "Invalid type for Scalar output.",
                        functionName,
                        __FILE__, __LINE__
                    );
                }

                if (sc1 < sc2)
                {
                    throw FoamXError
                    (
                        E_UNEXPECTED,
                        "Label value " + name(scalar(sc1))
                      + " is greater than the maximum " + name(scalar(sc2)),
                        functionName,
                        __FILE__, __LINE__
                    );
                }
            }
            break;

            default:
                throw FoamXError
                (
                    E_UNEXPECTED,
                    "Invalid type for operation <.  "
                    "Expected scalars or labels, found "
                  + FoamXTypes::typeName(any1.type),
                    functionName,
                    __FILE__, __LINE__
                );
        }
    }
    CATCH_ALL(functionName);
}


// ************************************************************************* //
