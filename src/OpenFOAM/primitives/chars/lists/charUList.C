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

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

static Ostream& writeChars
(
    Ostream& os,
    const char* chars,
    std::streamsize count
)
{
    // Contents are binary and contiguous
    os << nl << label(count) << nl;

    if (count)
    {
        const auto oldFmt = os.format(IOstream::BINARY);

        // write(...) includes surrounding start/end delimiters
        os.write(chars, count);

        os.format(oldFmt);
    }

    os.check(FUNCTION_NAME);
    return os;
}

} // End namespace Foam


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

namespace Foam
{

template<>
void UList<char>::writeEntry(Ostream& os) const
{
    const std::streamsize count(this->size());

    os  << word("List<char>");

    if (count)
    {
        writeChars(os, this->cdata(), count);
    }
    else
    {
        // Zero-sized binary - Write size only
        os << token::SPACE << label(0);
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace Foam
{

template<>
Ostream& UList<char>::writeList
(
    Ostream& os,
    const label /* shortLen (unused) */
) const
{
    return writeChars(os, this->cdata(), std::streamsize(this->size()));
}


template<>
Istream& UList<char>::readList(Istream& is)
{
    UList<char>& list = *this;

    // The target list length - must match with sizes read
    const label len = list.size();

    is.fatalCheck(FUNCTION_NAME);

    token tok(is);

    is.fatalCheck("UList<char>::readList(Istream&) : reading first token");

    if (tok.isCompound())
    {
        // Compound: simply transfer contents

        List<char> elems;
        elems.transfer
        (
            dynamicCast<token::Compound<List<char>>>
            (
                tok.transferCompoundToken(is)
            )
        );

        const label inputLen = elems.size();

        // List lengths must match
        if (inputLen != len)
        {
            FatalIOErrorInFunction(is)
                << "incorrect length for UList. Read "
                << inputLen << " expected " << len
                << exit(FatalIOError);
        }

        this->deepCopy(elems);
    }
    if (tok.isLabel())
    {
        // Label: could be int(..) or just a plain '0'

        const label inputLen = tok.labelToken();

        // List lengths must match
        if (inputLen != len)
        {
            FatalIOErrorInFunction(is)
                << "incorrect length for UList. Read "
                << inputLen << " expected " << len
                << exit(FatalIOError);
        }

        // Binary, always contiguous

        if (len)
        {
            const auto oldFmt = is.format(IOstream::BINARY);

            // read(...) includes surrounding start/end delimiters
            is.read(list.data(), std::streamsize(len));

            is.format(oldFmt);

            is.fatalCheck
            (
                "UList<char>::readList(Istream&) : "
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
