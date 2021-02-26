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

#include "charList.H"
#include "Istream.H"
#include "Ostream.H"
#include "token.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineCompoundTypeName(List<char>, charList);
    addCompoundToRunTimeSelectionTable(List<char>, charList);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace Foam
{

template<>
Istream& List<char>::readList(Istream& is)
{
    List<char>& list = *this;

    // Anull list
    list.clear();

    is.fatalCheck(FUNCTION_NAME);

    token tok(is);

    is.fatalCheck("List<char>::readList(Istream&) : reading first token");

    if (tok.isCompound())
    {
        // Compound: simply transfer contents

        list.transfer
        (
            dynamicCast<token::Compound<List<char>>>
            (
                tok.transferCompoundToken(is)
            )
        );
    }
    else if (tok.isLabel())
    {
        // Label: could be int(..) or just a plain '0'

        const label len = tok.labelToken();

        // Resize to actual length read
        list.resize(len);

        // Binary, always contiguous

        if (len)
        {
            const auto oldFmt = is.format(IOstream::BINARY);

            // read(...) includes surrounding start/end delimiters
            is.read(list.data(), std::streamsize(len));

            is.format(oldFmt);

            is.fatalCheck
            (
                "List<char>::readList(Istream&) : "
                "reading binary block"
            );
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "incorrect first token, expected <int>, found "
            << tok.info() << nl
            << exit(FatalIOError);
    }

    return is;
}


} // End namespace Foam


// ************************************************************************* //
