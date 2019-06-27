/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2019 OpenCFD Ltd.
     \\/     M anipulation  |
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, class Addr>
Foam::label Foam::IndirectListBase<T, Addr>::find
(
    const T& val,
    const label start
) const
{
    const label len = addr_.size();

    if (start >= 0 && len)
    {
        List_CONST_ACCESS(T, values_, vals);

        for (label i = start; i < len; ++i)
        {
            if (vals[addr_[i]] == val)
            {
                return i;
            }
        }
    }

    return -1;
}


template<class T, class Addr>
Foam::label Foam::IndirectListBase<T, Addr>::rfind
(
    const T& val,
    const label pos
) const
{
    List_CONST_ACCESS(T, values_, vals);

    const label len1 = (addr_.size()-1);

    // pos == -1 has same meaning as std::string::npos - search from end
    for (label i = ((pos >= 0 && pos < len1) ? pos : len1); i >= 0; --i)
    {
        if (vals[addr_[i]] == val)
        {
            return i;
        }
    }

    return -1;
}


// ************************************************************************* //
