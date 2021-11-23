/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "exprTraits.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::expressions::valueTypeCode
Foam::expressions::valueTypeCodeOf(const word& dataTypeName)
{
    #undef  stringToTypeCode
    #define stringToTypeCode(Type)                                  \
                                                                    \
    if (dataTypeName == exprTypeTraits<Type>::name)                 \
    {                                                               \
        return expressions::valueTypeCode::type_##Type;             \
    }

    if (!dataTypeName.empty())
    {
        stringToTypeCode(bool);
        stringToTypeCode(label);
        stringToTypeCode(scalar);
        stringToTypeCode(vector);
        stringToTypeCode(tensor);
        stringToTypeCode(sphericalTensor);
        stringToTypeCode(symmTensor);
    }
    #undef stringToTypeCode

    return expressions::valueTypeCode::INVALID;
}


Foam::word Foam::name(const expressions::valueTypeCode typeCode)
{
    #undef  case_typeCodeToString
    #define case_typeCodeToString(Type)                             \
                                                                    \
    case expressions::valueTypeCode::type_##Type :                  \
    {                                                               \
        return exprTypeTraits<Type>::name;                          \
    }

    switch (typeCode)
    {
        case expressions::valueTypeCode::NONE :
        {
            return "none";
        }

        case expressions::valueTypeCode::INVALID :
        {
            // ie, ""
            break;
        }

        case_typeCodeToString(bool);
        case_typeCodeToString(label);
        case_typeCodeToString(scalar);
        case_typeCodeToString(vector);
        case_typeCodeToString(tensor);
        case_typeCodeToString(sphericalTensor);
        case_typeCodeToString(symmTensor);
    }
    #undef case_typeCodeToString

    return word();
}


// ************************************************************************* //
