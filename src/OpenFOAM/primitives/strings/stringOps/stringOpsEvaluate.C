/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "stringOpsEvaluate.H"
#include "stringOps.H"
#include "StringStream.H"
#include "fieldExprDriver.H"
#include "error.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::string Foam::stringOps::evaluate
(
    const std::string& str,
    size_t pos,
    size_t len
)
{
    /// InfoErr<< "Evaluate " << str.substr(pos, len) << nl;

    size_t end = str.length();
    if (pos > end)
    {
        pos = end;
    }
    else if (len != std::string::npos)
    {
        len += pos;

        if (len < end)
        {
            end = len;
        }
    }

    // Adjust like inplaceTrim

    // Right
    while (pos < end && std::isspace(str[end-1]))
    {
        --end;
    }

    // Left
    while (pos < end && std::isspace(str[pos]))
    {
        ++pos;
    }

    if ((pos >= end) || std::isspace(str[pos]))
    {
        return "";
    }

    len = (end - pos);

    /// InfoErr<< "Evaluate " << str.substr(pos, len) << nl;

    expressions::exprResult result;
    {
        expressions::fieldExprDriver driver(1);
        driver.parse(str, pos, len);
        result = std::move(driver.result());
    }

    if (!result.hasValue() || !result.size())
    {
        InfoErr
            << "Failed evaluation: "
            << str.substr(pos, len) << nl;

        return "";
    }

    OStringStream os;
    result.writeValue(os);

    return os.str();
}


// ************************************************************************* //
