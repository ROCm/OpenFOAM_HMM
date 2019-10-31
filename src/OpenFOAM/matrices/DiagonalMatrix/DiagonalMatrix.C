/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "DiagonalMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
inline Foam::DiagonalMatrix<Type>::DiagonalMatrix()
:
    List<Type>()
{}


template<class Type>
Foam::DiagonalMatrix<Type>::DiagonalMatrix(const label n)
:
    List<Type>(n)
{}


template<class Type>
Foam::DiagonalMatrix<Type>::DiagonalMatrix(const label n, const zero)
:
    List<Type>(n, Zero)
{}


template<class Type>
Foam::DiagonalMatrix<Type>::DiagonalMatrix(const label n, const Type& val)
:
    List<Type>(n, val)
{}


template<class Type>
template<class Form>
Foam::DiagonalMatrix<Type>::DiagonalMatrix(const Matrix<Form, Type>& mat)
:
    List<Type>(min(mat.m(), mat.n()))
{
    label i = 0;

    for (Type& val : *this)
    {
        val = mat(i, i);
        ++i;
    }
}


template<class Type>
Foam::DiagonalMatrix<Type>& Foam::DiagonalMatrix<Type>::invert()
{
    for (Type& val : *this)
    {
        if (mag(val) < VSMALL)
        {
            val = Zero;
        }
        else
        {
            val = Type(1)/val;
        }
    }

    return this;
}


template<class Type>
Foam::DiagonalMatrix<Type> Foam::inv(const DiagonalMatrix<Type>& mat)
{
    DiagonalMatrix<Type> Ainv(mat.size());

    Type* iter = Ainv.begin();

    for (const Type& val : mat)
    {
        if (mag(val) < VSMALL)
        {
            *iter = Zero;
        }
        else
        {
            *iter = Type(1)/val;
        }

        ++iter;
    }

    return Ainv;
}


// ************************************************************************* //
