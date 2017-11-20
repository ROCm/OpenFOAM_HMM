/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2107 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
inline void Foam::vtk::write
(
    vtk::formatter& fmt,
    const Type& val
)
{
    const direction nCmpt = pTraits<Type>::nComponents;
    for (direction cmpt=0; cmpt < nCmpt; ++cmpt)
    {
        fmt.write(component(val, cmpt));
    }
}


template<class Type>
void Foam::vtk::writeList
(
    vtk::formatter& fmt,
    const UList<Type>& lst
)
{
    forAll(lst, i)
    {
        write(fmt, lst[i]);
    }
}


template<class Type, unsigned Size>
void Foam::vtk::writeList
(
    vtk::formatter& fmt,
    const FixedList<Type, Size>& lst
)
{
    for (unsigned i=0; i<Size; ++i)
    {
        write(fmt, lst[i]);
    }
}


template<class Type>
void Foam::vtk::writeList
(
    vtk::formatter& fmt,
    const UList<Type>& lst,
    const labelUList& addressing
)
{
    forAll(addressing, i)
    {
        write(fmt, lst[addressing[i]]);
    }
}


// ************************************************************************* //
