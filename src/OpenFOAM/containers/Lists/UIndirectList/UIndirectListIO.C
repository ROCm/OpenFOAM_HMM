/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

#include "UIndirectList.H"
#include "Ostream.H"
#include "token.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * Ostream Operator *  * * * * * * * * * * * * //

template<class T>
Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::UIndirectList<T>& L
)
{
    // Write list contents depending on data format
    if (os.format() == IOstream::ASCII || !contiguous<T>())
    {
        // Can the contents be considered 'uniform' (ie, identical)?
        bool uniform = (L.size() > 1 && contiguous<T>());
        if (uniform)
        {
            forAll(L, i)
            {
                if (L[i] != L[0])
                {
                    uniform = false;
                    break;
                }
            }
        }

        if (uniform)
        {
            // Write size and start delimiter
            os << L.size() << token::BEGIN_BLOCK;

            // Write contents
            os << L[0];

            // Write end delimiter
            os << token::END_BLOCK;
        }
        else if (L.size() <= 1 || (L.size() < 11 && contiguous<T>()))
        {
            // Write size and start delimiter
            os << L.size() << token::BEGIN_LIST;

            // Write contents
            forAll(L, i)
            {
                if (i) os << token::SPACE;
                os << L[i];
            }

            // Write end delimiter
            os << token::END_LIST;
        }
        else
        {
            // Write size and start delimiter
            os << nl << L.size() << nl << token::BEGIN_LIST;

            // Write contents
            forAll(L, i)
            {
                os << nl << L[i];
            }

            // Write end delimiter
            os << nl << token::END_LIST << nl;
        }
    }
    else
    {
        // Contents are binary and contiguous
        os << nl << L.size() << nl;

        if (L.size())
        {
            // This is annoying, and wasteful, but currently no alternative
            List<T> lst = L();

            // write(...) includes surrounding start/end delimiters
            os.write
            (
                reinterpret_cast<const char*>(lst.cdata()),
                lst.byteSize()
            );
        }
    }

    // Check state of IOstream
    os.check("Ostream& operator<<(Ostream&, const UIndirectList&)");

    return os;
}


// ************************************************************************* //
