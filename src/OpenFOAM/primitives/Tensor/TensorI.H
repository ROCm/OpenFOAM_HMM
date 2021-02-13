/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
#include <type_traits>

#include "SymmTensor.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Cmpt>
inline Foam::Tensor<Cmpt>::Tensor(const Foam::zero)
:
    Tensor::msType(Zero)
{}


template<class Cmpt>
template<class Cmpt2>
inline Foam::Tensor<Cmpt>::Tensor
(
    const MatrixSpace<Tensor<Cmpt2>, Cmpt2, 3, 3>& vs
)
:
    Tensor::msType(vs)
{}


template<class Cmpt>
template<class Cmpt2>
inline Foam::Tensor<Cmpt>::Tensor
(
    const VectorSpace<Tensor<Cmpt2>, Cmpt2, 9>& vs
)
:
    Tensor::msType(vs)
{}


template<class Cmpt>
inline Foam::Tensor<Cmpt>::Tensor(const SphericalTensor<Cmpt>& st)
{
    this->v_[XX] = st.ii();    this->v_[XY] = Zero;      this->v_[XZ] = Zero;
    this->v_[YX] = Zero;       this->v_[YY] = st.ii();   this->v_[YZ] = Zero;
    this->v_[ZX] = Zero;       this->v_[ZY] = Zero;      this->v_[ZZ] = st.ii();
}


template<class Cmpt>
inline Foam::Tensor<Cmpt>::Tensor(const SymmTensor<Cmpt>& st)
{
    this->v_[XX] = st.xx(); this->v_[XY] = st.xy(); this->v_[XZ] = st.xz();
    this->v_[YX] = st.xy(); this->v_[YY] = st.yy(); this->v_[YZ] = st.yz();
    this->v_[ZX] = st.xz(); this->v_[ZY] = st.yz(); this->v_[ZZ] = st.zz();
}


template<class Cmpt>
inline Foam::Tensor<Cmpt>::Tensor
(
    const Vector<Vector<Cmpt>>& vecs,
    const bool transposed
)
:
    Tensor<Cmpt>(vecs.x(), vecs.y(), vecs.z(), transposed)
{}


template<class Cmpt>
inline Foam::Tensor<Cmpt>::Tensor
(
    const Vector<Cmpt>& x,
    const Vector<Cmpt>& y,
    const Vector<Cmpt>& z,
    const bool transposed
)
{
    if (transposed)
    {
        this->cols(x, y, z);
    }
    else
    {
        this->rows(x, y, z);
    }
}


template<class Cmpt>
inline Foam::Tensor<Cmpt>::Tensor
(
    const Cmpt txx, const Cmpt txy, const Cmpt txz,
    const Cmpt tyx, const Cmpt tyy, const Cmpt tyz,
    const Cmpt tzx, const Cmpt tzy, const Cmpt tzz
)
{
    this->v_[XX] = txx; this->v_[XY] = txy; this->v_[XZ] = txz;
    this->v_[YX] = tyx; this->v_[YY] = tyy; this->v_[YZ] = tyz;
    this->v_[ZX] = tzx; this->v_[ZY] = tzy; this->v_[ZZ] = tzz;
}


template<class Cmpt>
template
<
    template<class, Foam::direction, Foam::direction> class Block2,
    Foam::direction BRowStart,
    Foam::direction BColStart
>
inline Foam::Tensor<Cmpt>::Tensor
(
    const Block2<Tensor<Cmpt>, BRowStart, BColStart>& block
)
:
    Tensor::msType(block)
{}


template<class Cmpt>
inline Foam::Tensor<Cmpt>::Tensor(Istream& is)
:
    Tensor::msType(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cmpt>
inline const Cmpt& Foam::Tensor<Cmpt>::xx() const
{
    return this->v_[XX];
}


template<class Cmpt>
inline const Cmpt& Foam::Tensor<Cmpt>::xy() const
{
    return this->v_[XY];
}


template<class Cmpt>
inline const Cmpt& Foam::Tensor<Cmpt>::xz() const
{
    return this->v_[XZ];
}


template<class Cmpt>
inline const Cmpt& Foam::Tensor<Cmpt>::yx() const
{
    return this->v_[YX];
}


template<class Cmpt>
inline const Cmpt& Foam::Tensor<Cmpt>::yy() const
{
    return this->v_[YY];
}


template<class Cmpt>
inline const Cmpt& Foam::Tensor<Cmpt>::yz() const
{
    return this->v_[YZ];
}


template<class Cmpt>
inline const Cmpt& Foam::Tensor<Cmpt>::zx() const
{
    return this->v_[ZX];
}


template<class Cmpt>
inline const Cmpt& Foam::Tensor<Cmpt>::zy() const
{
    return this->v_[ZY];
}


template<class Cmpt>
inline const Cmpt& Foam::Tensor<Cmpt>::zz() const
{
    return this->v_[ZZ];
}


template<class Cmpt>
inline Cmpt& Foam::Tensor<Cmpt>::xx()
{
    return this->v_[XX];
}


template<class Cmpt>
inline Cmpt& Foam::Tensor<Cmpt>::xy()
{
    return this->v_[XY];
}


template<class Cmpt>
inline Cmpt& Foam::Tensor<Cmpt>::xz()
{
    return this->v_[XZ];
}


template<class Cmpt>
inline Cmpt& Foam::Tensor<Cmpt>::yx()
{
    return this->v_[YX];
}


template<class Cmpt>
inline Cmpt& Foam::Tensor<Cmpt>::yy()
{
    return this->v_[YY];
}


template<class Cmpt>
inline Cmpt& Foam::Tensor<Cmpt>::yz()
{
    return this->v_[YZ];
}


template<class Cmpt>
inline Cmpt& Foam::Tensor<Cmpt>::zx()
{
    return this->v_[ZX];
}


template<class Cmpt>
inline Cmpt& Foam::Tensor<Cmpt>::zy()
{
    return this->v_[ZY];
}


template<class Cmpt>
inline Cmpt& Foam::Tensor<Cmpt>::zz()
{
    return this->v_[ZZ];
}


template<class Cmpt>
inline Foam::Vector<Cmpt> Foam::Tensor<Cmpt>::x() const
{
    return Vector<Cmpt>(this->v_[XX], this->v_[XY], this->v_[XZ]);
}


template<class Cmpt>
inline Foam::Vector<Cmpt> Foam::Tensor<Cmpt>::y() const
{
    return Vector<Cmpt>(this->v_[YX], this->v_[YY], this->v_[YZ]);
}


template<class Cmpt>
inline Foam::Vector<Cmpt> Foam::Tensor<Cmpt>::z() const
{
    return Vector<Cmpt>(this->v_[ZX], this->v_[ZY], this->v_[ZZ]);
}


template<class Cmpt>
inline Foam::Vector<Cmpt> Foam::Tensor<Cmpt>::cx() const
{
    return Vector<Cmpt>(this->v_[XX], this->v_[YX], this->v_[ZX]);
}


template<class Cmpt>
inline Foam::Vector<Cmpt> Foam::Tensor<Cmpt>::cy() const
{
    return Vector<Cmpt>(this->v_[XY], this->v_[YY], this->v_[ZY]);
}


template<class Cmpt>
inline Foam::Vector<Cmpt> Foam::Tensor<Cmpt>::cz() const
{
    return Vector<Cmpt>(this->v_[XZ], this->v_[YZ], this->v_[ZZ]);
}


template<class Cmpt>
template<Foam::direction Col>
inline Foam::Vector<Cmpt> Foam::Tensor<Cmpt>::col() const
{
    if (Col == 0) return cx();
    else if (Col == 1) return cy();
    else if (Col == 2) return cz();

    static_assert(Col < 3, "Invalid column access");
    return Zero;
}


template<class Cmpt>
inline Foam::Vector<Cmpt> Foam::Tensor<Cmpt>::col(const direction c) const
{
    switch (c)
    {
        case 0: return cx(); break;
        case 1: return cy(); break;
        case 2: return cz(); break;
    }

    FatalErrorInFunction
        << "Invalid column access " << c << abort(FatalError);

    return Zero;
}


template<class Cmpt>
template<Foam::direction Row>
inline Foam::Vector<Cmpt> Foam::Tensor<Cmpt>::row() const
{
    if (Row == 0) return x();
    else if (Row == 1) return y();
    else if (Row == 2) return z();

    static_assert(Row < 3, "Invalid row access");
    return Zero;
}


template<class Cmpt>
inline Foam::Vector<Cmpt> Foam::Tensor<Cmpt>::row(const direction r) const
{
    switch (r)
    {
        case 0: return x(); break;
        case 1: return y(); break;
        case 2: return z(); break;
    }

    FatalErrorInFunction
        << "Invalid row access " << r << abort(FatalError);

    return Zero;
}


template<class Cmpt>
template<Foam::direction Col>
inline void Foam::Tensor<Cmpt>::col(const Vector<Cmpt>& v)
{
    if (Col == 0)
    {
        this->v_[XX] = v.x();
        this->v_[YX] = v.y();
        this->v_[ZX] = v.z();
    }
    else if (Col == 1)
    {
        this->v_[XY] = v.x();
        this->v_[YY] = v.y();
        this->v_[ZY] = v.z();
    }
    else if (Col == 2)
    {
        this->v_[XZ] = v.x();
        this->v_[YZ] = v.y();
        this->v_[ZZ] = v.z();
    }

    static_assert(Col < 3, "Invalid column access");
}


template<class Cmpt>
template<Foam::direction Row>
inline void Foam::Tensor<Cmpt>::row(const Vector<Cmpt>& v)
{
    if (Row == 0)
    {
        this->v_[XX] = v.x(); this->v_[XY] = v.y(); this->v_[XZ] = v.z();
    }
    else if (Row == 1)
    {
        this->v_[YX] = v.x(); this->v_[YY] = v.y(); this->v_[YZ] = v.z();
    }
    else if (Row == 2)
    {
        this->v_[ZX] = v.x(); this->v_[ZY] = v.y(); this->v_[ZZ] = v.z();
    }

    static_assert(Row < 3, "Invalid row access");
}


template<class Cmpt>
inline void Foam::Tensor<Cmpt>::cols
(
    const Vector<Cmpt>& x,
    const Vector<Cmpt>& y,
    const Vector<Cmpt>& z
)
{
    this->v_[XX] = x.x(); this->v_[XY] = y.x(); this->v_[XZ] = z.x();
    this->v_[YX] = x.y(); this->v_[YY] = y.y(); this->v_[YZ] = z.y();
    this->v_[ZX] = x.z(); this->v_[ZY] = y.z(); this->v_[ZZ] = z.z();
}


template<class Cmpt>
inline void Foam::Tensor<Cmpt>::rows
(
    const Vector<Cmpt>& x,
    const Vector<Cmpt>& y,
    const Vector<Cmpt>& z
)
{
    this->v_[XX] = x.x(); this->v_[XY] = x.y(); this->v_[XZ] = x.z();
    this->v_[YX] = y.x(); this->v_[YY] = y.y(); this->v_[YZ] = y.z();
    this->v_[ZX] = z.x(); this->v_[ZY] = z.y(); this->v_[ZZ] = z.z();
}


template<class Cmpt>
inline void Foam::Tensor<Cmpt>::col(const direction c, const Vector<Cmpt>& v)
{
    switch (c)
    {
        case 0: col<0>(v); break;
        case 1: col<1>(v); break;
        case 2: col<2>(v); break;
        default:
            FatalErrorInFunction
                << "Invalid column access " << c << abort(FatalError);
            break;
    }
}


template<class Cmpt>
inline void Foam::Tensor<Cmpt>::row(const direction r, const Vector<Cmpt>& v)
{
    switch (r)
    {
        case 0: row<0>(v); break;
        case 1: row<1>(v); break;
        case 2: row<2>(v); break;
        default:
            FatalErrorInFunction
                << "Invalid row access " << r << abort(FatalError);
            break;
    }
}


template<class Cmpt>
inline Foam::Vector<Cmpt> Foam::Tensor<Cmpt>::diag() const
{
    return Vector<Cmpt>(this->v_[XX], this->v_[YY], this->v_[ZZ]);
}


template<class Cmpt>
inline void Foam::Tensor<Cmpt>::diag(const Vector<Cmpt>& v)
{
    this->v_[XX] = v.x(); this->v_[YY] = v.y(); this->v_[ZZ] = v.z();
}


// * * * * * * * * * * * * * * * Member Operations * * * * * * * * * * * * * //

template<class Cmpt>
inline Foam::Tensor<Cmpt> Foam::Tensor<Cmpt>::T() const
{
    return Tensor<Cmpt>
    (
        xx(), yx(), zx(),
        xy(), yy(), zy(),
        xz(), yz(), zz()
    );
}


template<class Cmpt>
inline Foam::Tensor<Cmpt>
Foam::Tensor<Cmpt>::inner(const Tensor<Cmpt>& t2) const
{
    const Tensor<Cmpt>& t1 = *this;

    return Tensor<Cmpt>
    (
        t1.xx()*t2.xx() + t1.xy()*t2.yx() + t1.xz()*t2.zx(),
        t1.xx()*t2.xy() + t1.xy()*t2.yy() + t1.xz()*t2.zy(),
        t1.xx()*t2.xz() + t1.xy()*t2.yz() + t1.xz()*t2.zz(),

        t1.yx()*t2.xx() + t1.yy()*t2.yx() + t1.yz()*t2.zx(),
        t1.yx()*t2.xy() + t1.yy()*t2.yy() + t1.yz()*t2.zy(),
        t1.yx()*t2.xz() + t1.yy()*t2.yz() + t1.yz()*t2.zz(),

        t1.zx()*t2.xx() + t1.zy()*t2.yx() + t1.zz()*t2.zx(),
        t1.zx()*t2.xy() + t1.zy()*t2.yy() + t1.zz()*t2.zy(),
        t1.zx()*t2.xz() + t1.zy()*t2.yz() + t1.zz()*t2.zz()
    );
}


template<class Cmpt>
inline Foam::Tensor<Cmpt>
Foam::Tensor<Cmpt>::schur(const Tensor<Cmpt>& t2) const
{
    const Tensor<Cmpt>& t1 = *this;

    return Tensor<Cmpt>
    (
        t1.xx()*t2.xx(), t1.xy()*t2.xy(), t1.xz()*t2.xz(),
        t1.yx()*t2.yx(), t1.yy()*t2.yy(), t1.yz()*t2.yz(),
        t1.zx()*t2.zx(), t1.zy()*t2.zy(), t1.zz()*t2.zz()
    );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Cmpt>
inline void Foam::Tensor<Cmpt>::operator&=(const Tensor<Cmpt>& t)
{
    *this = this->inner(t);
}


template<class Cmpt>
template<class Cmpt2>
inline void Foam::Tensor<Cmpt>::operator=
(
    const VectorSpace<Tensor<Cmpt2>, Cmpt2, 9>& vs
)
{
    VectorSpace<Tensor<Cmpt>, Cmpt, 9>::operator=(vs);
}


template<class Cmpt>
inline void Foam::Tensor<Cmpt>::operator=(const SphericalTensor<Cmpt>& st)
{
    this->v_[XX] = st.ii();  this->v_[XY] = Zero;     this->v_[XZ] = Zero;
    this->v_[YX] = Zero;     this->v_[YY] = st.ii();  this->v_[YZ] = Zero;
    this->v_[ZX] = Zero;     this->v_[ZY] = Zero;     this->v_[ZZ] = st.ii();
}


template<class Cmpt>
inline void Foam::Tensor<Cmpt>::operator=(const SymmTensor<Cmpt>& st)
{
    this->v_[XX] = st.xx(); this->v_[XY] = st.xy(); this->v_[XZ] = st.xz();
    this->v_[YX] = st.xy(); this->v_[YY] = st.yy(); this->v_[YZ] = st.yz();
    this->v_[ZX] = st.xz(); this->v_[ZY] = st.yz(); this->v_[ZZ] = st.zz();
}


template<class Cmpt>
inline void Foam::Tensor<Cmpt>::operator=(const Vector<Vector<Cmpt>>& tr)
{
    this->v_[XX] = tr.x().x();
    this->v_[XY] = tr.x().y();
    this->v_[XZ] = tr.x().z();

    this->v_[YX] = tr.y().x();
    this->v_[YY] = tr.y().y();
    this->v_[YZ] = tr.y().z();

    this->v_[ZX] = tr.z().x();
    this->v_[ZY] = tr.z().y();
    this->v_[ZZ] = tr.z().z();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Return the trace of a Tensor
template<class Cmpt>
inline Cmpt tr(const Tensor<Cmpt>& t)
{
    return t.xx() + t.yy() + t.zz();
}


//- Return the spherical part of a Tensor
template<class Cmpt>
inline SphericalTensor<Cmpt> sph(const Tensor<Cmpt>& t)
{
    return SphericalTensor<Cmpt>
    (
        (1.0/3.0)*tr(t)
    );
}


//- Return the symmetric part of a Tensor
template<class Cmpt>
inline SymmTensor<Cmpt> symm(const Tensor<Cmpt>& t)
{
    return SymmTensor<Cmpt>
    (
        t.xx(), 0.5*(t.xy() + t.yx()), 0.5*(t.xz() + t.zx()),
                t.yy(),                0.5*(t.yz() + t.zy()),
                                       t.zz()
    );
}


//- Return twice the symmetric part of a Tensor
template<class Cmpt>
inline SymmTensor<Cmpt> twoSymm(const Tensor<Cmpt>& t)
{
    return SymmTensor<Cmpt>
    (
        2*t.xx(), (t.xy() + t.yx()), (t.xz() + t.zx()),
                  2*t.yy(),          (t.yz() + t.zy()),
                                     2*t.zz()
    );
}


//- Return the skew-symmetric part of a Tensor
template<class Cmpt>
inline Tensor<Cmpt> skew(const Tensor<Cmpt>& t)
{
    return Tensor<Cmpt>
    (
        Zero,                  0.5*(t.xy() - t.yx()), 0.5*(t.xz() - t.zx()),
        0.5*(t.yx() - t.xy()), Zero,                  0.5*(t.yz() - t.zy()),
        0.5*(t.zx() - t.xz()), 0.5*(t.zy() - t.yz()), Zero
    );
}


//- Return the skew-symmetric part of a SymmTensor as a Tensor
template<class Cmpt>
inline const Tensor<Cmpt>& skew(const SymmTensor<Cmpt>& st)
{
    return Tensor<Cmpt>::zero;
}


//- Return the deviatoric part of a Tensor
template<class Cmpt>
inline Tensor<Cmpt> dev(const Tensor<Cmpt>& t)
{
    return t - sph(t);
}


//- Return the two-third deviatoric part of a Tensor
template<class Cmpt>
inline Tensor<Cmpt> dev2(const Tensor<Cmpt>& t)
{
    return t - 2*sph(t);
}


//- Return the determinant of a Tensor
template<class Cmpt>
inline Cmpt det(const Tensor<Cmpt>& t)
{
    return
    (
        t.xx()*t.yy()*t.zz() + t.xy()*t.yz()*t.zx()
      + t.xz()*t.yx()*t.zy() - t.xx()*t.yz()*t.zy()
      - t.xy()*t.yx()*t.zz() - t.xz()*t.yy()*t.zx()
    );
}


//- Return the cofactor Tensor of a Tensor
template<class Cmpt>
inline Tensor<Cmpt> cof(const Tensor<Cmpt>& t)
{
    return Tensor<Cmpt>
    (
        t.yy()*t.zz() - t.zy()*t.yz(),
        t.zx()*t.yz() - t.yx()*t.zz(),
        t.yx()*t.zy() - t.yy()*t.zx(),

        t.xz()*t.zy() - t.xy()*t.zz(),
        t.xx()*t.zz() - t.xz()*t.zx(),
        t.xy()*t.zx() - t.xx()*t.zy(),

        t.xy()*t.yz() - t.xz()*t.yy(),
        t.yx()*t.xz() - t.xx()*t.yz(),
        t.xx()*t.yy() - t.yx()*t.xy()
    );
}


//- Return the inverse of a Tensor by using the given determinant
template<class Cmpt>
inline Tensor<Cmpt> inv(const Tensor<Cmpt>& t, const Cmpt dett)
{
    #ifdef FULLDEBUG
    if (mag(dett) < VSMALL)
    {
        FatalErrorInFunction
            << "Tensor is not invertible due to the (almost) zero determinant:"
            << "    Tensor = " << t << nl
            << "    det(Tensor) = " << dett
            << abort(FatalError);
    }
    #endif

    return cof(t).T()/dett;
}


//- Return the inverse of a Tensor
template<class Cmpt>
inline Tensor<Cmpt> inv(const Tensor<Cmpt>& t)
{
    return inv(t, det(t));
}


//- Return the inverse of this Tensor
template<class Cmpt>
inline Tensor<Cmpt> Tensor<Cmpt>::inv() const
{
    return Foam::inv(*this);
}


//- Return the 1st invariant of a Tensor
template<class Cmpt>
inline Cmpt invariantI(const Tensor<Cmpt>& t)
{
    return tr(t);
}


//- Return the 2nd invariant of a Tensor
template<class Cmpt>
inline Cmpt invariantII(const Tensor<Cmpt>& t)
{
    return
    (
        t.xx()*t.yy() + t.yy()*t.zz() + t.xx()*t.zz()
      - t.xy()*t.yx() - t.yz()*t.zy() - t.xz()*t.zx()
    );
}


//- Return the 3rd invariant of a Tensor
template<class Cmpt>
inline Cmpt invariantIII(const Tensor<Cmpt>& t)
{
    return det(t);
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Sum of a SphericalTensor and a Tensor
template<class Cmpt>
inline Tensor<Cmpt>
operator+(const SphericalTensor<Cmpt>& st1, const Tensor<Cmpt>& t2)
{
    return Tensor<Cmpt>
    (
        st1.ii() + t2.xx(), t2.xy(),            t2.xz(),
        t2.yx(),            st1.ii() + t2.yy(), t2.yz(),
        t2.zx(),            t2.zy(),            st1.ii() + t2.zz()
    );
}


//- Sum of a Tensor and a SphericalTensor
template<class Cmpt>
inline Tensor<Cmpt>
operator+(const Tensor<Cmpt>& t1, const SphericalTensor<Cmpt>& st2)
{
    return Tensor<Cmpt>
    (
        t1.xx() + st2.ii(), t1.xy(),            t1.xz(),
        t1.yx(),            t1.yy() + st2.ii(), t1.yz(),
        t1.zx(),            t1.zy(),            t1.zz() + st2.ii()
    );
}


//- Sum of a SymmTensor and a Tensor
template<class Cmpt>
inline Tensor<Cmpt>
operator+(const SymmTensor<Cmpt>& st1, const Tensor<Cmpt>& t2)
{
    return Tensor<Cmpt>
    (
        st1.xx() + t2.xx(), st1.xy() + t2.xy(), st1.xz() + t2.xz(),
        st1.xy() + t2.yx(), st1.yy() + t2.yy(), st1.yz() + t2.yz(),
        st1.xz() + t2.zx(), st1.yz() + t2.zy(), st1.zz() + t2.zz()
    );
}


//- Sum of a Tensor and a SymmTensor
template<class Cmpt>
inline Tensor<Cmpt>
operator+(const Tensor<Cmpt>& t1, const SymmTensor<Cmpt>& st2)
{
    return Tensor<Cmpt>
    (
        t1.xx() + st2.xx(), t1.xy() + st2.xy(), t1.xz() + st2.xz(),
        t1.yx() + st2.xy(), t1.yy() + st2.yy(), t1.yz() + st2.yz(),
        t1.zx() + st2.xz(), t1.zy() + st2.yz(), t1.zz() + st2.zz()
    );
}


//- Subtract a Tensor from a SphericalTensor
template<class Cmpt>
inline Tensor<Cmpt>
operator-(const SphericalTensor<Cmpt>& st1, const Tensor<Cmpt>& t2)
{
    return Tensor<Cmpt>
    (
        st1.ii() - t2.xx(), -t2.xy(),            -t2.xz(),
       -t2.yx(),             st1.ii() - t2.yy(), -t2.yz(),
       -t2.zx(),            -t2.zy(),             st1.ii() - t2.zz()
    );
}


//- Subtract a SphericalTensor from a Tensor
template<class Cmpt>
inline Tensor<Cmpt>
operator-(const Tensor<Cmpt>& t1, const SphericalTensor<Cmpt>& st2)
{
    return Tensor<Cmpt>
    (
        t1.xx() - st2.ii(), t1.xy(),            t1.xz(),
        t1.yx(),            t1.yy() - st2.ii(), t1.yz(),
        t1.zx(),            t1.zy(),            t1.zz() - st2.ii()
    );
}


//- Subtract a Tensor from a SymmTensor
template<class Cmpt>
inline Tensor<Cmpt>
operator-(const SymmTensor<Cmpt>& st1, const Tensor<Cmpt>& t2)
{
    return Tensor<Cmpt>
    (
        st1.xx() - t2.xx(), st1.xy() - t2.xy(), st1.xz() - t2.xz(),
        st1.xy() - t2.yx(), st1.yy() - t2.yy(), st1.yz() - t2.yz(),
        st1.xz() - t2.zx(), st1.yz() - t2.zy(), st1.zz() - t2.zz()
    );
}


//- Subtract a SymmTensor from a Tensor
template<class Cmpt>
inline Tensor<Cmpt>
operator-(const Tensor<Cmpt>& t1, const SymmTensor<Cmpt>& st2)
{
    return Tensor<Cmpt>
    (
        t1.xx() - st2.xx(), t1.xy() - st2.xy(), t1.xz() - st2.xz(),
        t1.yx() - st2.xy(), t1.yy() - st2.yy(), t1.yz() - st2.yz(),
        t1.zx() - st2.xz(), t1.zy() - st2.yz(), t1.zz() - st2.zz()
    );
}


//- Return the Hodge dual of a Tensor as a Vector
template<class Cmpt>
inline Vector<Cmpt> operator*(const Tensor<Cmpt>& t)
{
    return Vector<Cmpt>(t.yz(), -t.xz(), t.xy());
}


//- Return the Hodge dual of a Vector as a Tensor
template<class Cmpt>
inline Tensor<Cmpt> operator*(const Vector<Cmpt>& v)
{
    return Tensor<Cmpt>
    (
         Zero,  -v.z(),  v.y(),
         v.z(),  Zero,  -v.x(),
        -v.y(),  v.x(),  Zero
    );
}


//- Division of a Vector by a Tensor
template<class Cmpt>
inline typename innerProduct<Vector<Cmpt>, Tensor<Cmpt>>::type
operator/(const Vector<Cmpt>& v, const Tensor<Cmpt>& t)
{
    return inv(t) & v;
}


//- Division of a Tensor by a Cmpt
template<class Cmpt>
inline Tensor<Cmpt>
operator/(const Tensor<Cmpt>& t, const Cmpt s)
{
    #ifdef FULLDEBUG
    if (mag(s) < VSMALL)
    {
        FatalErrorInFunction
            << "Tensor = " << t
            << " is not divisible due to a zero value in Cmpt:"
            << "Cmpt = " << s
            << abort(FatalError);
    }
    #endif

    return Tensor<Cmpt>
    (
        t.xx()/s, t.xy()/s, t.xz()/s,
        t.yx()/s, t.yy()/s, t.yz()/s,
        t.zx()/s, t.zy()/s, t.zz()/s
    );
}


//- Inner-product of a Tensor and a Tensor
template<class Cmpt>
inline typename innerProduct<Tensor<Cmpt>, Tensor<Cmpt>>::type
operator&(const Tensor<Cmpt>& t1, const Tensor<Cmpt>& t2)
{
    return t1.inner(t2);
}


//- Inner-product of a SphericalTensor and a Tensor
template<class Cmpt>
inline Tensor<Cmpt>
operator&(const SphericalTensor<Cmpt>& st1, const Tensor<Cmpt>& t2)
{
    return Tensor<Cmpt>
    (
        st1.ii()*t2.xx(), st1.ii()*t2.xy(), st1.ii()*t2.xz(),
        st1.ii()*t2.yx(), st1.ii()*t2.yy(), st1.ii()*t2.yz(),
        st1.ii()*t2.zx(), st1.ii()*t2.zy(), st1.ii()*t2.zz()
    );
}


//- Inner-product of a Tensor and a SphericalTensor
template<class Cmpt>
inline Tensor<Cmpt>
operator&(const Tensor<Cmpt>& t1, const SphericalTensor<Cmpt>& st2)
{
    return Tensor<Cmpt>
    (
        t1.xx()*st2.ii(), t1.xy()*st2.ii(), t1.xz()*st2.ii(),
        t1.yx()*st2.ii(), t1.yy()*st2.ii(), t1.yz()*st2.ii(),
        t1.zx()*st2.ii(), t1.zy()*st2.ii(), t1.zz()*st2.ii()
    );
}


//- Inner-product of a SymmTensor and a Tensor
template<class Cmpt>
inline Tensor<Cmpt>
operator&(const SymmTensor<Cmpt>& st1, const Tensor<Cmpt>& t2)
{
    return Tensor<Cmpt>
    (
        st1.xx()*t2.xx() + st1.xy()*t2.yx() + st1.xz()*t2.zx(),
        st1.xx()*t2.xy() + st1.xy()*t2.yy() + st1.xz()*t2.zy(),
        st1.xx()*t2.xz() + st1.xy()*t2.yz() + st1.xz()*t2.zz(),

        st1.xy()*t2.xx() + st1.yy()*t2.yx() + st1.yz()*t2.zx(),
        st1.xy()*t2.xy() + st1.yy()*t2.yy() + st1.yz()*t2.zy(),
        st1.xy()*t2.xz() + st1.yy()*t2.yz() + st1.yz()*t2.zz(),

        st1.xz()*t2.xx() + st1.yz()*t2.yx() + st1.zz()*t2.zx(),
        st1.xz()*t2.xy() + st1.yz()*t2.yy() + st1.zz()*t2.zy(),
        st1.xz()*t2.xz() + st1.yz()*t2.yz() + st1.zz()*t2.zz()
    );
}


//- Inner-product of a Tensor and a SymmTensor
template<class Cmpt>
inline Tensor<Cmpt>
operator&(const Tensor<Cmpt>& t1, const SymmTensor<Cmpt>& st2)
{
    return Tensor<Cmpt>
    (
        t1.xx()*st2.xx() + t1.xy()*st2.xy() + t1.xz()*st2.xz(),
        t1.xx()*st2.xy() + t1.xy()*st2.yy() + t1.xz()*st2.yz(),
        t1.xx()*st2.xz() + t1.xy()*st2.yz() + t1.xz()*st2.zz(),

        t1.yx()*st2.xx() + t1.yy()*st2.xy() + t1.yz()*st2.xz(),
        t1.yx()*st2.xy() + t1.yy()*st2.yy() + t1.yz()*st2.yz(),
        t1.yx()*st2.xz() + t1.yy()*st2.yz() + t1.yz()*st2.zz(),

        t1.zx()*st2.xx() + t1.zy()*st2.xy() + t1.zz()*st2.xz(),
        t1.zx()*st2.xy() + t1.zy()*st2.yy() + t1.zz()*st2.yz(),
        t1.zx()*st2.xz() + t1.zy()*st2.yz() + t1.zz()*st2.zz()
    );
}


//- Inner-product of a Tensor and a Vector
template<class Cmpt>
inline typename innerProduct<Tensor<Cmpt>, Vector<Cmpt>>::type
operator&(const Tensor<Cmpt>& t, const Vector<Cmpt>& v)
{
    return Vector<Cmpt>
    (
        t.xx()*v.x() + t.xy()*v.y() + t.xz()*v.z(),
        t.yx()*v.x() + t.yy()*v.y() + t.yz()*v.z(),
        t.zx()*v.x() + t.zy()*v.y() + t.zz()*v.z()
    );
}


//- Inner-product of a Vector and a Tensor
template<class Cmpt>
inline typename innerProduct<Vector<Cmpt>, Tensor<Cmpt>>::type
operator&(const Vector<Cmpt>& v, const Tensor<Cmpt>& t)
{
    return Vector<Cmpt>
    (
        v.x()*t.xx() + v.y()*t.yx() + v.z()*t.zx(),
        v.x()*t.xy() + v.y()*t.yy() + v.z()*t.zy(),
        v.x()*t.xz() + v.y()*t.yz() + v.z()*t.zz()
    );
}


//- Double-inner-product of a SphericalTensor and a Tensor
template<class Cmpt>
inline Cmpt
operator&&(const SphericalTensor<Cmpt>& st1, const Tensor<Cmpt>& t2)
{
    return (st1.ii()*t2.xx() + st1.ii()*t2.yy() + st1.ii()*t2.zz());
}


//- Double-inner-product of a Tensor and a SphericalTensor
template<class Cmpt>
inline Cmpt
operator&&(const Tensor<Cmpt>& t1, const SphericalTensor<Cmpt>& st2)
{
    return (t1.xx()*st2.ii() + t1.yy()*st2.ii() + t1.zz()*st2.ii());
}


//- Double-inner-product of a SymmTensor and a Tensor
template<class Cmpt>
inline Cmpt
operator&&(const SymmTensor<Cmpt>& st1, const Tensor<Cmpt>& t2)
{
    return
    (
        st1.xx()*t2.xx() + st1.xy()*t2.xy() + st1.xz()*t2.xz() +
        st1.xy()*t2.yx() + st1.yy()*t2.yy() + st1.yz()*t2.yz() +
        st1.xz()*t2.zx() + st1.yz()*t2.zy() + st1.zz()*t2.zz()
    );
}


//- Double-inner-product of a Tensor and a SymmTensor
template<class Cmpt>
inline Cmpt
operator&&(const Tensor<Cmpt>& t1, const SymmTensor<Cmpt>& st2)
{
    return
    (
        t1.xx()*st2.xx() + t1.xy()*st2.xy() + t1.xz()*st2.xz() +
        t1.yx()*st2.xy() + t1.yy()*st2.yy() + t1.yz()*st2.yz() +
        t1.zx()*st2.xz() + t1.zy()*st2.yz() + t1.zz()*st2.zz()
    );
}


//- Outer-product of a Vector and a Vector
template<class Cmpt>
inline typename outerProduct<Vector<Cmpt>, Vector<Cmpt>>::type
operator*(const Vector<Cmpt>& v1, const Vector<Cmpt>& v2)
{
    return Tensor<Cmpt>
    (
        v1.x()*v2.x(), v1.x()*v2.y(), v1.x()*v2.z(),
        v1.y()*v2.x(), v1.y()*v2.y(), v1.y()*v2.z(),
        v1.z()*v2.x(), v1.z()*v2.y(), v1.z()*v2.z()
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Cmpt>
class typeOfSum<SphericalTensor<Cmpt>, Tensor<Cmpt>>
{
public:

    typedef Tensor<Cmpt> type;
};


template<class Cmpt>
class typeOfSum<Tensor<Cmpt>, SphericalTensor<Cmpt>>
{
public:

    typedef Tensor<Cmpt> type;
};


template<class Cmpt>
class innerProduct<SphericalTensor<Cmpt>, Tensor<Cmpt>>
{
public:

    typedef Tensor<Cmpt> type;
};


template<class Cmpt>
class innerProduct<Tensor<Cmpt>, SphericalTensor<Cmpt>>
{
public:

    typedef Tensor<Cmpt> type;
};


template<class Cmpt>
class typeOfSum<SymmTensor<Cmpt>, Tensor<Cmpt>>
{
public:

    typedef Tensor<Cmpt> type;
};


template<class Cmpt>
class typeOfSum<Tensor<Cmpt>, SymmTensor<Cmpt>>
{
public:

    typedef Tensor<Cmpt> type;
};


template<class Cmpt>
class innerProduct<SymmTensor<Cmpt>, Tensor<Cmpt>>
{
public:

    typedef Tensor<Cmpt> type;
};


template<class Cmpt>
class innerProduct<Tensor<Cmpt>, SymmTensor<Cmpt>>
{
public:

    typedef Tensor<Cmpt> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
