/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

template<class T>
inline Foam::label Foam::UIndirectList<T>::find
(
    const T& val,
    const label start
) const
{
    if (start >= 0)
    {
        List_CONST_ACCESS(T, completeList_, lst);
        List_CONST_ACCESS(label, addressing_, addr);

        const label len = addressing_.size();

        for (label i = start; i < len; ++i)
        {
            if (lst[addr[i]] == val)
            {
                return i;
            }
        }
    }

    return -1;
}


template<class T>
inline Foam::label Foam::UIndirectList<T>::rfind
(
    const T& val,
    const label pos
) const
{
    List_CONST_ACCESS(T, completeList_, lst);
    List_CONST_ACCESS(label, addressing_, addr);

    for
    (
        label i =
        (
            pos < 0
          ? (addressing_.size()-1)
          : min(pos, (addressing_.size()-1))
        );
        i >= 0;
        --i
    )
    {
        if (lst[addr[i]] == val)
        {
            return i;
        }
    }

    return -1;
}


// ************************************************************************* //
