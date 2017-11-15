/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class T>
Foam::Ostream& Foam::UIndirectList<T>::writeList
(
    Ostream& os,
    const label shortListLen
) const
{
    const UIndirectList<T>& L = *this;

    const label sz = L.size();

    // Write list contents depending on data format
    if (os.format() == IOstream::ASCII || !contiguous<T>())
    {
        // Can the contents be considered 'uniform' (ie, identical)?
        bool uniform = (sz > 1 && contiguous<T>());
        if (uniform)
        {
            for (label i=1; i < sz; ++i)
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
            // Size and start delimiter
            os << sz << token::BEGIN_BLOCK;

            // Contents
            os << L[0];

            // End delimiter
            os << token::END_BLOCK;
        }
        else if
        (
            sz <= 1 || !shortListLen
         || (sz <= shortListLen && contiguous<T>())
        )
        {
            // Size and start delimiter
            os << sz << token::BEGIN_LIST;

            // Contents
            for (label i=0; i < sz; ++i)
            {
                if (i) os << token::SPACE;
                os << L[i];
            }

            // End delimiter
            os << token::END_LIST;
        }
        else
        {
            // Size and start delimiter
            os << nl << sz << nl << token::BEGIN_LIST << nl;

            // Contents
            for (label i=0; i < sz; ++i)
            {
                os << L[i] << nl;
            }

            // End delimiter
            os << token::END_LIST << nl;
        }
    }
    else
    {
        // Contents are binary and contiguous
        os << nl << sz << nl;

        if (sz)
        {
            // The TOTAL number of bytes to be written.
            // - possibly add start delimiter
            os.beginRaw(sz*sizeof(T));

            // Contents
            for (label i=0; i < sz; ++i)
            {
                os.writeRaw
                (
                    reinterpret_cast<const char*>(&(L[i])),
                    sizeof(T)
                );
            }

            // End delimiter and/or cleanup.
            os.endRaw();
        }
    }

    os.check(FUNCTION_NAME);
    return os;
}


// * * * * * * * * * * * * * * * Ostream Operator *  * * * * * * * * * * * * //

template<class T>
Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::UIndirectList<T>& lst
)
{
    return lst.writeList(os, 10);
}


// ************************************************************************* //
