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
    const UIndirectList<T>& list = *this;

    const label len = list.size();

    // Write list contents depending on data format
    if (os.format() == IOstream::ASCII || !contiguous<T>())
    {
        if (contiguous<T>() && list.uniform())
        {
            // Two or more entries, and all entries have identical values.
            os << len << token::BEGIN_BLOCK << list[0] << token::END_BLOCK;
        }
        else if
        (
            len <= 1 || !shortListLen
         || (len <= shortListLen && contiguous<T>())
        )
        {
            // Size and start delimiter
            os << len << token::BEGIN_LIST;

            // Contents
            for (label i=0; i < len; ++i)
            {
                if (i) os << token::SPACE;
                os << list[i];
            }

            // End delimiter
            os << token::END_LIST;
        }
        else
        {
            // Size and start delimiter
            os << nl << len << nl << token::BEGIN_LIST << nl;

            // Contents
            for (label i=0; i < len; ++i)
            {
                os << list[i] << nl;
            }

            // End delimiter
            os << token::END_LIST << nl;
        }
    }
    else
    {
        // Contents are binary and contiguous
        os << nl << len << nl;

        if (len)
        {
            // The TOTAL number of bytes to be written.
            // - possibly add start delimiter
            os.beginRaw(len*sizeof(T));

            // Contents
            for (label i=0; i < len; ++i)
            {
                os.writeRaw
                (
                    reinterpret_cast<const char*>(&(list[i])),
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
    const Foam::UIndirectList<T>& list
)
{
    return list.writeList(os, 10);
}


// ************************************************************************* //
