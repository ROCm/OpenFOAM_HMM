/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "MatrixBlock.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Form, class Type>
inline void Foam::Matrix<Form, Type>::doAlloc()
{
    const label len = size();

    if (len)
    {
        v_ = new Type[len];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Form, class Type>
inline Foam::Matrix<Form, Type>::Matrix()
:
    mRows_(0),
    nCols_(0),
    v_(nullptr)
{}


template<class Form, class Type>
inline Foam::autoPtr<Foam::Matrix<Form, Type>>
Foam::Matrix<Form, Type>::clone() const
{
    return autoPtr<Matrix<Form, Type>>::New(*this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Form, class Type>
inline const Foam::Matrix<Form, Type>& Foam::Matrix<Form, Type>::null()
{
    return NullObjectRef<Matrix<Form, Type>>();
}


template<class Form, class Type>
inline Foam::label Foam::Matrix<Form, Type>::m() const noexcept
{
    return mRows_;
}


template<class Form, class Type>
inline Foam::label Foam::Matrix<Form, Type>::n() const noexcept
{
    return nCols_;
}


template<class Form, class Type>
inline Foam::label Foam::Matrix<Form, Type>::size() const
{
    return mRows_ * nCols_;
}


template<class Form, class Type>
inline bool Foam::Matrix<Form, Type>::empty() const noexcept
{
    return !mRows_ || !nCols_;
}


template<class Form, class Type>
inline void Foam::Matrix<Form, Type>::checki(const label i) const
{
    #ifdef FULLDEBUG
    if (!mRows_ || !nCols_)
    {
        FatalErrorInFunction
            << "Attempt to access element from empty matrix"
            << abort(FatalError);
    }
    if (i < 0 || i >= mRows_)
    {
        FatalErrorInFunction
            << "Index " << i << " out of range 0 ... " << mRows_-1
            << abort(FatalError);
    }
    #endif
}


template<class Form, class Type>
inline void Foam::Matrix<Form, Type>::checkj(const label j) const
{
    #ifdef FULLDEBUG
    if (!mRows_ || !nCols_)
    {
        FatalErrorInFunction
            << "Attempt to access element from empty matrix"
            << abort(FatalError);
    }
    if (j < 0 || j >= nCols_)
    {
        FatalErrorInFunction
            << "index " << j << " out of range 0 ... " << nCols_-1
            << abort(FatalError);
    }
    #endif
}


template<class Form, class Type>
inline void Foam::Matrix<Form, Type>::checkSize() const
{
    if (mRows_ < 0 || nCols_ < 0)
    {
        FatalErrorInFunction
            << "Incorrect size (" << mRows_ << ", " << nCols_ << ')' << nl
            << abort(FatalError);
    }
    // Could also check for odd sizes, like (0, N) and make (0,0)
}


template<class Form, class Type>
inline bool Foam::Matrix<Form, Type>::uniform() const
{
    const label len = size();

    if (len == 0)
    {
        return false;
    }

    for (label idx=1; idx<len; ++idx)
    {
        if (v_[0] != v_[idx])
        {
            return false;
        }
    }

    return true;
}


template<class Form, class Type>
inline const Type* Foam::Matrix<Form, Type>::cdata() const noexcept
{
    return v_;
}


template<class Form, class Type>
inline Type* Foam::Matrix<Form, Type>::data() noexcept
{
    return v_;
}


template<class Form, class Type>
inline const Type* Foam::Matrix<Form, Type>::rowData(const label irow) const
{
    #ifdef FULLDEBUG
    checki(irow);
    #endif
    return (v_ + irow*nCols_);
}


template<class Form, class Type>
inline Type* Foam::Matrix<Form, Type>::rowData(const label irow)
{
    #ifdef FULLDEBUG
    checki(irow);
    #endif
    return (v_ + irow*nCols_);
}


template<class Form, class Type>
inline const Type& Foam::Matrix<Form, Type>::at(const label idx) const
{
    #ifdef FULLDEBUG
    if (idx < 0 || idx >= this->size())
    {
        FatalErrorInFunction
            << "index " << idx << " out of range 0 ... " << this->size()
            << abort(FatalError);
    }
    #endif
    return (v_ + idx);
}


template<class Form, class Type>
inline Type& Foam::Matrix<Form, Type>::at(const label idx)
{
    #ifdef FULLDEBUG
    if (idx < 0 || idx >= this->size())
    {
        FatalErrorInFunction
            << "index " << idx << " out of range 0 ... " << this->size()
            << abort(FatalError);
    }
    #endif
    return (v_ + idx);
}


template<class Form, class Type>
inline Foam::ConstMatrixBlock<Foam::Matrix<Form, Type>>
Foam::Matrix<Form, Type>::block
(
    const label m,
    const label n,
    const label mStart,
    const label nStart
) const
{
    return ConstMatrixBlock<mType>
    (
        *this,
        m,
        n,
        mStart,
        nStart
    );
}


template<class Form, class Type>
template<class VectorSpace>
inline Foam::ConstMatrixBlock<Foam::Matrix<Form, Type>>
Foam::Matrix<Form, Type>::block
(
    const label mStart,
    const label nStart
) const
{
    return ConstMatrixBlock<mType>
    (
        *this,
        VectorSpace::mRows,
        VectorSpace::nCols,
        mStart,
        nStart
    );
}


template<class Form, class Type>
inline Foam::ConstMatrixBlock<Foam::Matrix<Form, Type>>
Foam::Matrix<Form, Type>::col
(
    const label m,
    const label mStart
) const
{
    return ConstMatrixBlock<mType>
    (
        *this,
        m,
        1,
        mStart,
        0
    );
}


template<class Form, class Type>
inline Foam::ConstMatrixBlock<Foam::Matrix<Form, Type>>
Foam::Matrix<Form, Type>::col
(
    const label m,
    const label mStart,
    const label nStart
) const
{
    return ConstMatrixBlock<mType>
    (
        *this,
        m,
        1,
        mStart,
        nStart
    );
}


template<class Form, class Type>
inline Foam::MatrixBlock<Foam::Matrix<Form, Type>>
Foam::Matrix<Form, Type>::block
(
    const label m,
    const label n,
    const label mStart,
    const label nStart
)
{
    return MatrixBlock<mType>
    (
        *this,
        m,
        n,
        mStart,
        nStart
    );
}


template<class Form, class Type>
template<class VectorSpace>
inline Foam::MatrixBlock<Foam::Matrix<Form, Type>>
Foam::Matrix<Form, Type>::block(const label mStart, const label nStart)
{
    return MatrixBlock<mType>
    (
        *this,
        VectorSpace::mRows,
        VectorSpace::nCols,
        mStart,
        nStart
    );
}


template<class Form, class Type>
inline Foam::MatrixBlock<Foam::Matrix<Form, Type>>
Foam::Matrix<Form, Type>::col(const label m, const label mStart)
{
    return MatrixBlock<mType>
    (
        *this,
        m,
        1,
        mStart,
        0
    );
}


template<class Form, class Type>
inline Foam::MatrixBlock<Foam::Matrix<Form, Type>>
Foam::Matrix<Form, Type>::col
(
    const label m,
    const label mStart,
    const label nStart
)
{
    return MatrixBlock<mType>
    (
        *this,
        m,
        1,
        mStart,
        nStart
    );
}


template<class Form, class Type>
inline void Foam::Matrix<Form, Type>::setSize(const label m, const label n)
{
    resize(m, n);
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::shallowResize(const label m, const label n)
{
    mRows_ = m;
    nCols_ = n;
}


template<class Form, class Type>
inline Foam::tmp<Foam::Field<Type>> Foam::Matrix<Form, Type>::rightMultiply
(
    const UList<Type>& colVec
) const
{
    return this->rightMultiplyImpl(colVec);
}


template<class Form, class Type>
template<class Addr>
inline Foam::tmp<Foam::Field<Type>> Foam::Matrix<Form, Type>::rightMultiply
(
    const IndirectListBase<Type, Addr>& colVec
) const
{
    return this->rightMultiplyImpl(colVec);
}


template<class Form, class Type>
inline Foam::tmp<Foam::Field<Type>> Foam::Matrix<Form, Type>::leftMultiply
(
    const UList<Type>& rowVec
) const
{
    return this->leftMultiplyImpl(rowVec);
}


template<class Form, class Type>
template<class Addr>
inline Foam::tmp<Foam::Field<Type>> Foam::Matrix<Form, Type>::leftMultiply
(
    const IndirectListBase<Type, Addr>& rowVec
) const
{
    return this->leftMultiplyImpl(rowVec);
}


// * * * * * * * * * * * * * * * * * Iterators * * * * * * * * * * * * * * * //

template<class Form, class Type>
inline typename Foam::Matrix<Form, Type>::iterator
Foam::Matrix<Form, Type>::begin()
{
    return v_;
}


template<class Form, class Type>
inline typename Foam::Matrix<Form, Type>::iterator
Foam::Matrix<Form, Type>::end()
{
    return v_ + (mRows_ * nCols_);
}


template<class Form, class Type>
inline typename Foam::Matrix<Form, Type>::const_iterator
Foam::Matrix<Form, Type>::cbegin() const
{
    return v_;
}


template<class Form, class Type>
inline typename Foam::Matrix<Form, Type>::const_iterator
Foam::Matrix<Form, Type>::cend() const
{
    return v_ + (mRows_ * nCols_);
}


template<class Form, class Type>
inline typename Foam::Matrix<Form, Type>::const_iterator
Foam::Matrix<Form, Type>::begin() const
{
    return v_;
}


template<class Form, class Type>
inline typename Foam::Matrix<Form, Type>::const_iterator
Foam::Matrix<Form, Type>::end() const
{
    return v_ + (mRows_ * nCols_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Form, class Type>
inline const Type& Foam::Matrix<Form, Type>::operator()
(
    const label irow,
    const label jcol
) const
{
    checki(irow);
    checkj(jcol);
    return v_[irow*nCols_ + jcol];
}


template<class Form, class Type>
inline Type& Foam::Matrix<Form, Type>::operator()
(
    const label irow,
    const label jcol
)
{
    checki(irow);
    checkj(jcol);
    return v_[irow*nCols_ + jcol];
}


template<class Form, class Type>
inline const Type* Foam::Matrix<Form, Type>::operator[](const label irow) const
{
    checki(irow);
    return v_ + irow*nCols_;
}


template<class Form, class Type>
inline Type* Foam::Matrix<Form, Type>::operator[](const label irow)
{
    checki(irow);
    return v_ + irow*nCols_;
}


template<class Form, class Type>
inline Foam::tmp<Foam::Field<Type>> Foam::operator*
(
    const Matrix<Form, Type>& mat,
    const UList<Type>& x
)
{
    mat.rightMultiply(x);
}


template<class Form, class Type, class Addr>
inline Foam::tmp<Foam::Field<Type>> Foam::operator*
(
    const Matrix<Form, Type>& mat,
    const IndirectListBase<Type, Addr>& x
)
{
    mat.rightMultiply(x);
}


template<class Form, class Type>
inline Foam::tmp<Foam::Field<Type>> Foam::operator*
(
    const UList<Type>& x,
    const Matrix<Form, Type>& mat
)
{
    mat.leftMultiply(x);
}


template<class Form, class Type, class Addr>
inline Foam::tmp<Foam::Field<Type>> Foam::operator*
(
    const IndirectListBase<Type, Addr>& x,
    const Matrix<Form, Type>& mat
)
{
    mat.leftMultiply(x);
}


// ************************************************************************* //
