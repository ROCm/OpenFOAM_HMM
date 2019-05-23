/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "Matrix.H"
#include <functional>
#include <algorithm>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Form, class Type>
template<class ListType>
Foam::tmp<Foam::Field<Type>> Foam::Matrix<Form, Type>::AmulImpl
(
    const ListType& x
) const
{
    const Matrix<Form, Type>& mat = *this;

    #ifdef FULLDEBUG
    if (mat.n() != x.size())
    {
        FatalErrorInFunction
            << "Attempt to multiply incompatible Matrix and Vector:" << nl
            << "Matrix : (" << mat.m() << ", " << mat.n() << ')' << nl
            << "Matrix columns != Vector size (" << x.size() << ')' << nl
            << abort(FatalError);
    }
    #endif

    auto tresult = tmp<Field<Type>>::New(mat.m(), Zero);
    auto& result = tresult.ref();

    for (label i = 0; i < mat.m(); ++i)
    {
        for (label j = 0; j < mat.n(); ++j)
        {
            result[i] += mat(i, j)*x[j];
        }
    }

    return tresult;
}


template<class Form, class Type>
template<class ListType>
Foam::tmp<Foam::Field<Type>> Foam::Matrix<Form, Type>::TmulImpl
(
    const ListType& x
) const
{
    const Matrix<Form, Type>& mat = *this;

    #ifdef FULLDEBUG
    if (mat.m() != x.size())
    {
        FatalErrorInFunction
            << "Attempt to multiply incompatible Matrix and Vector:" << nl
            << "Matrix : (" << mat.m() << ", " << mat.n() << ')' << nl
            << "Matrix rows != Vector size (" << x.size() << ')' << nl
            << abort(FatalError);
    }
    #endif

    auto tresult = tmp<Field<Type>>::New(mat.n(), Zero);
    auto& result = tresult.ref();

    for (label i = 0; i < mat.m(); ++i)
    {
        const Type& val = x[i];
        for (label j = 0; j < mat.n(); ++j)
        {
            result[j] += val*mat(i, j);
        }
    }

    return tresult;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix(const label m, const label n)
:
    mRows_(m),
    nCols_(n),
    v_(nullptr)
{
    checkSize();

    doAlloc();
}


template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix(const label m, const label n, const zero)
:
    mRows_(m),
    nCols_(n),
    v_(nullptr)
{
    checkSize();

    doAlloc();

    std::fill(begin(), end(), Zero);
}


template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix(const label m, const label n, const Type& val)
:
    mRows_(m),
    nCols_(n),
    v_(nullptr)
{
    checkSize();

    doAlloc();

    std::fill(begin(), end(), val);
}


template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix(const Matrix<Form, Type>& mat)
:
    mRows_(mat.mRows_),
    nCols_(mat.nCols_),
    v_(nullptr)
{
    if (mat.cdata())
    {
        doAlloc();

        std::copy(mat.cbegin(), mat.cend(), v_);
    }
}


template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix(Matrix<Form, Type>&& mat)
:
    mRows_(mat.mRows_),
    nCols_(mat.nCols_),
    v_(mat.v_)
{
    mat.mRows_ = 0;
    mat.nCols_ = 0;
    mat.v_ = nullptr;
}


template<class Form, class Type>
template<class Form2>
Foam::Matrix<Form, Type>::Matrix(const Matrix<Form2, Type>& mat)
:
    mRows_(mat.m()),
    nCols_(mat.n()),
    v_(nullptr)
{
    if (mat.cdata())
    {
        doAlloc();

        std::copy(mat.cbegin(), mat.cend(), v_);
    }
}


template<class Form, class Type>
template<class MatrixType>
inline Foam::Matrix<Form, Type>::Matrix
(
    const ConstMatrixBlock<MatrixType>& Mb
)
:
    mRows_(Mb.m()),
    nCols_(Mb.n())
{
    doAlloc();

    for (label i = 0; i < mRows_; ++i)
    {
        for (label j = 0; j < nCols_; ++j)
        {
            (*this)(i, j) = Mb(i,j);
        }
    }
}


template<class Form, class Type>
template<class MatrixType>
inline Foam::Matrix<Form, Type>::Matrix
(
    const MatrixBlock<MatrixType>& Mb
)
:
    mRows_(Mb.m()),
    nCols_(Mb.n())
{
    doAlloc();

    for (label i = 0; i < mRows_; ++i)
    {
        for (label j = 0; j < nCols_; ++j)
        {
            (*this)(i, j) = Mb(i, j);
        }
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class Form, class Type>
Foam::Matrix<Form, Type>::~Matrix()
{
    if (v_)
    {
        delete[] v_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Form, class Type>
void Foam::Matrix<Form, Type>::clear()
{
    if (v_)
    {
        delete[] v_;
        v_ = nullptr;
    }

    mRows_ = 0;
    nCols_ = 0;
}


template<class Form, class Type>
Foam::List<Type> Foam::Matrix<Form, Type>::release()
{
    List<Type> list;

    const label len = size();

    if (v_ && len)
    {
        UList<Type> storage(v_, len);
        list.swap(storage);

        v_ = nullptr;
    }
    clear();

    return list;
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::swap(Matrix<Form, Type>& mat)
{
    Foam::Swap(mRows_, mat.mRows_);
    Foam::Swap(nCols_, mat.nCols_);
    Foam::Swap(v_, mat.v_);
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::transfer(Matrix<Form, Type>& mat)
{
    clear();

    mRows_ = mat.mRows_;
    nCols_ = mat.nCols_;
    v_ = mat.v_;

    mat.mRows_ = 0;
    mat.nCols_ = 0;
    mat.v_ = nullptr;
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::resize(const label m, const label n)
{
    if (m == mRows_ && n == nCols_)
    {
        return;
    }

    Matrix<Form, Type> newMatrix(m, n, Zero);

    const label mrow = min(m, mRows_);
    const label ncol = min(n, nCols_);

    for (label i = 0; i < mrow; ++i)
    {
        for (label j = 0; j < ncol; ++j)
        {
            newMatrix(i, j) = (*this)(i, j);
        }
    }

    transfer(newMatrix);
}


template<class Form, class Type>
Form Foam::Matrix<Form, Type>::T() const
{
    Form At(n(), m());

    for (label i = 0; i < m(); ++i)
    {
        for (label j = 0; j < n(); ++j)
        {
            At(j, i) = (*this)(i, j);
        }
    }

    return At;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator=(const Matrix<Form, Type>& mat)
{
    if (this == &mat)
    {
        FatalErrorInFunction
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    if (mRows_ != mat.mRows_ || nCols_ != mat.nCols_)
    {
        clear();
        mRows_ = mat.mRows_;
        nCols_ = mat.nCols_;
        doAlloc();
    }

    if (v_)
    {
        std::copy(mat.cbegin(), mat.cend(), v_);
    }
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator=(Matrix<Form, Type>&& mat)
{
    if (this == &mat)
    {
        FatalErrorInFunction
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    this->transfer(mat);
}


template<class Form, class Type>
template<class MatrixType>
void Foam::Matrix<Form, Type>::operator=
(
    const ConstMatrixBlock<MatrixType>& Mb
)
{
    for (label i = 0; i < mRows_; ++i)
    {
        for (label j = 0; j < nCols_; ++j)
        {
            (*this)(i, j) = Mb(i, j);
        }
    }
}


template<class Form, class Type>
template<class MatrixType>
void Foam::Matrix<Form, Type>::operator=
(
    const MatrixBlock<MatrixType>& Mb
)
{
    for (label i=0; i < mRows_; ++i)
    {
        for (label j=0; j < nCols_; ++j)
        {
            (*this)(i, j) = Mb(i, j);
        }
    }
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator=(const Type& val)
{
    std::fill(begin(), end(), val);
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator=(const zero)
{
    std::fill(begin(), end(), Zero);
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator+=(const Matrix<Form, Type>& other)
{
    if (this == &other)
    {
        FatalErrorInFunction
            << "Attempted addition to self"
            << abort(FatalError);
    }

    if (m() != other.m() || n() != other.n())
    {
        FatalErrorInFunction
            << "Attempt to add matrices with different sizes: ("
            << m() << ", " << n() << ") != ("
            << other.m() << ", " << other.n() << ')' << nl
            << abort(FatalError);
    }

    Type* out = this->data();
    const Type* in = other.cdata();

    const label len = this->size();

    for (label idx = 0; idx < len; ++idx)
    {
        out[idx] += in[idx];
    }
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator-=(const Matrix<Form, Type>& other)
{
    if (this == &other)
    {
        FatalErrorInFunction
            << "Attempted subtraction from self"
            << abort(FatalError);
    }

    if (m() != other.m() || n() != other.n())
    {
        FatalErrorInFunction
            << "Attempt to subtract matrices with different sizes: ("
            << m() << ", " << n() << ") != ("
            << other.m() << ", " << other.n() << ')' << nl
            << abort(FatalError);
    }

    Type* out = this->data();
    const Type* in = other.cdata();

    const label len = this->size();

    for (label idx=0; idx < len; ++idx)
    {
        out[idx] -= in[idx];
    }
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator*=(const scalar s)
{
    for (Type& val : *this)
    {
        val *= s;
    }
}


template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator/=(const scalar s)
{
    for (Type& val : *this)
    {
        val /= s;
    }
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Form, class Type>
const Type& Foam::max(const Matrix<Form, Type>& mat)
{
    if (mat.empty())
    {
        FatalErrorInFunction
            << "Matrix is empty" << abort(FatalError);
    }

    return *(std::max_element(mat.cbegin(), mat.cend()));
}


template<class Form, class Type>
const Type& Foam::min(const Matrix<Form, Type>& mat)
{
    if (mat.empty())
    {
        FatalErrorInFunction
            << "Matrix is empty" << abort(FatalError);
    }

    return *(std::min_element(mat.cbegin(), mat.cend()));
}


template<class Form, class Type>
Foam::MinMax<Type> Foam::minMax(const Matrix<Form, Type>& mat)
{
    MinMax<Type> result;

    for (const Type& val : mat)
    {
        result += val;
    }

    return result;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Form, class Type>
Form Foam::operator-(const Matrix<Form, Type>& mat)
{
    Form result(mat.m(), mat.n());

    std::transform
    (
        mat.cbegin(),
        mat.cend(),
        result.begin(),
        std::negate<Type>()
    );

    return result;
}


template<class Form, class Type>
Form Foam::operator+(const Matrix<Form, Type>& A, const Matrix<Form, Type>& B)
{
    if (A.m() != B.m())
    {
        FatalErrorInFunction
            << "Attempt to add matrices with different sizes: ("
            << A.m() << ", " << A.n() << ") != ("
            << B.m() << ", " << B.n() << ')' << nl
            << abort(FatalError);
    }

    Form AB(A.m(), A.n());

    Type* ABv = AB.data();
    const Type* Av = A.cdata();
    const Type* Bv = B.cdata();

    const label len = A.size();

    for (label idx = 0; idx < len; ++idx)
    {
        ABv[idx] = Av[idx] + Bv[idx];
    }

    return AB;
}


template<class Form, class Type>
Form Foam::operator-(const Matrix<Form, Type>& A, const Matrix<Form, Type>& B)
{
    if (A.m() != B.m())
    {
        FatalErrorInFunction
            << "Attempt to subtract matrices with different sizes: ("
            << A.m() << ", " << A.n() << ") != ("
            << B.m() << ", " << B.n() << ')' << nl
            << abort(FatalError);
    }

    Form AB(A.m(), A.n());

    Type* ABv = AB.data();
    const Type* Av = A.cdata();
    const Type* Bv = B.cdata();

    const label len = A.size();

    for (label idx=0; idx < len; ++idx)
    {
        ABv[idx] = Av[idx] - Bv[idx];
    }

    return AB;
}


template<class Form, class Type>
Form Foam::operator*(const scalar s, const Matrix<Form, Type>& mat)
{
    Form result(mat.m(), mat.n());

    const label len = mat.size();

    if (len)
    {
        Type* out = result.data();
        const Type* in = mat.cdata();

        for (label idx = 0; idx < len; ++idx)
        {
            out[idx] = s * in[idx];
        }
    }

    return result;
}


template<class Form, class Type>
Form Foam::operator*(const Matrix<Form, Type>& mat, const scalar s)
{
    Form result(mat.m(), mat.n());

    const label len = mat.size();

    if (len)
    {
        Type* out = result.data();
        const Type* in = mat.cdata();

        for (label idx=0; idx < len; ++idx)
        {
            out[idx] = in[idx] * s;
        }
    }

    return result;
}


template<class Form, class Type>
Form Foam::operator/(const Matrix<Form, Type>& mat, const scalar s)
{
    Form result(mat.m(), mat.n());

    const label len = mat.size();

    if (len)
    {
        Type* out = result.data();
        const Type* in = mat.cdata();

        for (label idx=0; idx < len; ++idx)
        {
            out[idx] = in[idx] / s;
        }
    }

    return result;
}


template<class Form1, class Form2, class Type>
typename Foam::typeOfInnerProduct<Type, Form1, Form2>::type
Foam::operator*
(
    const Matrix<Form1, Type>& A,
    const Matrix<Form2, Type>& B
)
{
    if (A.n() != B.m())
    {
        FatalErrorInFunction
            << "Attempt to multiply incompatible matrices:" << nl
            << "Matrix A : (" << A.m() << ", " << A.n() << ')' << nl
            << "Matrix B : (" << B.m() << ", " << B.n() << ')' << nl
            << "The columns of A must equal rows of B"
            << abort(FatalError);
    }

    typename typeOfInnerProduct<Type, Form1, Form2>::type AB
    (
        A.m(),
        B.n(),
        Zero
    );

    for (label i=0; i < AB.m(); ++i)
    {
        for (label j=0; j < AB.n(); ++j)
        {
            for (label k=0; k < B.m(); ++k)
            {
                AB(i, j) += A(i, k) * B(k, j);
            }
        }
    }

    return AB;
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "MatrixIO.C"

// ************************************************************************* //
