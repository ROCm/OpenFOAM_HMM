/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
inline Foam::Matrix<Type>::Matrix()
:
    n_(0),
    m_(0),
    v_(NULL)
{}


template<class Type>
inline Foam::autoPtr<Foam::Matrix<Type> > Foam::Matrix<Type>::clone() const
{
    return autoPtr<Matrix<Type> >(new Matrix<Type>(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return the number of rows
template<class Type>
inline Foam::label Foam::Matrix<Type>::n() const
{
    return n_;
}


template<class Type>
inline Foam::label Foam::Matrix<Type>::m() const
{
    return m_;
}


template<class Type>
inline Foam::label Foam::Matrix<Type>::size() const
{
    return n_*m_;
}


template<class Type>
inline void Foam::Matrix<Type>::checki(const label i) const
{
    if (!n_)
    {
        FatalErrorIn("Matrix<Type>::checki(const label)")
            << "attempt to access element from zero sized row"
            << abort(FatalError);
    }
    else if (i<0 || i>=n_)
    {
        FatalErrorIn("Matrix<Type>::checki(const label)")
            << "index " << i << " out of range 0 ... " << n_-1
            << abort(FatalError);
    }
}


template<class Type>
inline void Foam::Matrix<Type>::checkj(const label j) const
{
    if (!m_)
    {
        FatalErrorIn("Matrix<Type>::checkj(const label)")
            << "attempt to access element from zero sized column"
            << abort(FatalError);
    }
    else if (j<0 || j>=m_)
    {
        FatalErrorIn("Matrix<Type>::checkj(const label)")
            << "index " << j << " out of range 0 ... " << m_-1
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
inline Type* Foam::Matrix<Type>::operator[](const label i)
{
#   ifdef FULLDEBUG
    checki(i);
#   endif
    return v_[i];
}


template<class Type>
inline const Type* Foam::Matrix<Type>::operator[](const label i) const
{
#   ifdef FULLDEBUG
    checki(i);
#   endif
    return v_[i];
}


// ************************************************************************* //
