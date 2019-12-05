/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "PtrListDetail.H"
#include "error.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
Foam::Ostream& Foam::Detail::PtrListDetail<T>::write
(
    Ostream& os,
    const bool trimNull
) const
{
    const label len = this->size();

    // The (output) size and start delimiter
    os  << nl << indent << (trimNull ? this->count() : len) << nl
        << indent << token::BEGIN_LIST << incrIndent << nl;

    // Contents
    for (label i=0; i < len; ++i)
    {
        const T* ptr = (*this)[i];
        if (ptr)
        {
            os << *ptr << nl;
        }
        else if (!trimNull)
        {
            FatalErrorInFunction
                << "cannot dereference nullptr at index " << i
                << " in range [0," << len << ")"
                << abort(FatalError);
        }
    }

    // End delimiter
    os << decrIndent << indent << token::END_LIST << nl;

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
