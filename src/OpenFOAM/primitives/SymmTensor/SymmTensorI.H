/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "Tensor.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Cmpt>
inline Foam::SymmTensor<Cmpt>::SymmTensor(const Foam::zero)
:
    SymmTensor::vsType(Zero)
{}


template<class Cmpt>
template<class Cmpt2>
inline Foam::SymmTensor<Cmpt>::SymmTensor
(
    const VectorSpace<SymmTensor<Cmpt2>, Cmpt2, 6>& vs
)
:
    SymmTensor::vsType(vs)
{}


template<class Cmpt>
inline Foam::SymmTensor<Cmpt>::SymmTensor(const SphericalTensor<Cmpt>& st)
{
    this->v_[XX] = st.ii(); this->v_[XY] = Zero; this->v_[XZ] = Zero;
    this->v_[YY] = st.ii(); this->v_[YZ] = Zero;
    this->v_[ZZ] = st.ii();
}


template<class Cmpt>
inline Foam::SymmTensor<Cmpt>::SymmTensor
(
    const Cmpt txx, const Cmpt txy, const Cmpt txz,
                    const Cmpt tyy, const Cmpt tyz,
                                    const Cmpt tzz
)
{
    this->v_[XX] = txx; this->v_[XY] = txy; this->v_[XZ] = txz;
                        this->v_[YY] = tyy; this->v_[YZ] = tyz;
                                            this->v_[ZZ] = tzz;
}


template<class Cmpt>
inline Foam::SymmTensor<Cmpt>::SymmTensor(Istream& is)
:
    SymmTensor::vsType(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cmpt>
inline const Cmpt& Foam::SymmTensor<Cmpt>::xx() const
{
    return this->v_[XX];
}

template<class Cmpt>
inline const Cmpt& Foam::SymmTensor<Cmpt>::xy() const
{
    return this->v_[XY];
}

template<class Cmpt>
inline const Cmpt& Foam::SymmTensor<Cmpt>::xz() const
{
    return this->v_[XZ];
}

template<class Cmpt>
inline const Cmpt& Foam::SymmTensor<Cmpt>::yx() const
{
    return this->v_[XY];
}

template<class Cmpt>
inline const Cmpt& Foam::SymmTensor<Cmpt>::yy() const
{
    return this->v_[YY];
}

template<class Cmpt>
inline const Cmpt& Foam::SymmTensor<Cmpt>::yz() const
{
    return this->v_[YZ];
}

template<class Cmpt>
inline const Cmpt& Foam::SymmTensor<Cmpt>::zx() const
{
    return this->v_[XZ];
}

template<class Cmpt>
inline const Cmpt& Foam::SymmTensor<Cmpt>::zy() const
{
    return this->v_[YZ];
}

template<class Cmpt>
inline const Cmpt& Foam::SymmTensor<Cmpt>::zz() const
{
    return this->v_[ZZ];
}


template<class Cmpt>
inline Cmpt& Foam::SymmTensor<Cmpt>::xx()
{
    return this->v_[XX];
}

template<class Cmpt>
inline Cmpt& Foam::SymmTensor<Cmpt>::xy()
{
    return this->v_[XY];
}

template<class Cmpt>
inline Cmpt& Foam::SymmTensor<Cmpt>::xz()
{
    return this->v_[XZ];
}

template<class Cmpt>
inline Cmpt& Foam::SymmTensor<Cmpt>::yx()
{
    return this->v_[XY];
}

template<class Cmpt>
inline Cmpt& Foam::SymmTensor<Cmpt>::yy()
{
    return this->v_[YY];
}

template<class Cmpt>
inline Cmpt& Foam::SymmTensor<Cmpt>::yz()
{
    return this->v_[YZ];
}

template<class Cmpt>
inline Cmpt& Foam::SymmTensor<Cmpt>::zx()
{
    return this->v_[XZ];
}

template<class Cmpt>
inline Cmpt& Foam::SymmTensor<Cmpt>::zy()
{
    return this->v_[YZ];
}

template<class Cmpt>
inline Cmpt& Foam::SymmTensor<Cmpt>::zz()
{
    return this->v_[ZZ];
}


template<class Cmpt>
inline Foam::Vector<Cmpt> Foam::SymmTensor<Cmpt>::diag() const
{
    return Vector<Cmpt>(this->v_[XX], this->v_[YY], this->v_[ZZ]);
}


template<class Cmpt>
inline void Foam::SymmTensor<Cmpt>::diag(const Vector<Cmpt>& v)
{
    this->v_[XX] = v.x(); this->v_[YY] = v.y(); this->v_[ZZ] = v.z();
}


template<class Cmpt>
inline const Foam::SymmTensor<Cmpt>& Foam::SymmTensor<Cmpt>::T() const
{
    return *this;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Cmpt>
inline void Foam::SymmTensor<Cmpt>::operator=(const SphericalTensor<Cmpt>& st)
{
    this->v_[XX] = st.ii(); this->v_[XY] = Zero;    this->v_[XZ] = Zero;
                            this->v_[YY] = st.ii(); this->v_[YZ] = Zero;
                                                    this->v_[ZZ] = st.ii();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Return the trace of a SymmTensor
template<class Cmpt>
inline Cmpt tr(const SymmTensor<Cmpt>& st)
{
    return st.xx() + st.yy() + st.zz();
}


//- Return the spherical part of a SymmTensor
template<class Cmpt>
inline SphericalTensor<Cmpt> sph(const SymmTensor<Cmpt>& st)
{
    return SphericalTensor<Cmpt>
    (
        (1.0/3.0)*tr(st)
    );
}


//- Return the symmetric part of a SymmTensor, i.e. itself
template<class Cmpt>
inline const SymmTensor<Cmpt>& symm(const SymmTensor<Cmpt>& st)
{
    return st;
}


//- Return twice the symmetric part of a SymmTensor
template<class Cmpt>
inline SymmTensor<Cmpt> twoSymm(const SymmTensor<Cmpt>& st)
{
    return 2*st;
}


//- Return the deviatoric part of a SymmTensor
template<class Cmpt>
inline SymmTensor<Cmpt> dev(const SymmTensor<Cmpt>& st)
{
    return st - sph(st);
}


//- Return the two-third deviatoric part of a SymmTensor
template<class Cmpt>
inline SymmTensor<Cmpt> dev2(const SymmTensor<Cmpt>& st)
{
    return st - 2*sph(st);
}


//- Return the determinant of a SymmTensor
template<class Cmpt>
inline Cmpt det(const SymmTensor<Cmpt>& st)
{
    return
    (
        st.xx()*st.yy()*st.zz() + st.xy()*st.yz()*st.xz()
      + st.xz()*st.xy()*st.yz() - st.xx()*st.yz()*st.yz()
      - st.xy()*st.xy()*st.zz() - st.xz()*st.yy()*st.xz()
    );
}


//- Return the cofactor of a SymmTensor as a SymmTensor
template<class Cmpt>
inline SymmTensor<Cmpt> cof(const SymmTensor<Cmpt>& st)
{
    return SymmTensor<Cmpt>
    (
        st.yy()*st.zz() - st.yz()*st.yz(),
        st.xz()*st.yz() - st.xy()*st.zz(),
        st.xy()*st.yz() - st.xz()*st.yy(),

        st.xx()*st.zz() - st.xz()*st.xz(),
        st.xy()*st.xz() - st.xx()*st.yz(),

        st.xx()*st.yy() - st.xy()*st.xy()
    );
}


//- Return the inverse of a SymmTensor by using the given determinant
template<class Cmpt>
inline SymmTensor<Cmpt> inv(const SymmTensor<Cmpt>& st, const Cmpt detst)
{
    #ifdef FULLDEBUG
    if (mag(detst) < VSMALL)
    {
        FatalErrorInFunction
            << "SymmTensor is not invertible due to the zero determinant:"
            << "det(symmTensor) = " << mag(detst)
            << abort(FatalError);
    }
    #endif

    return cof(st).T()/detst;
}


//- Return the inverse of a SymmTensor
template<class Cmpt>
inline SymmTensor<Cmpt> inv(const SymmTensor<Cmpt>& st)
{
    return inv(st, det(st));
}


//- Return the 1st invariant of a SymmTensor
template<class Cmpt>
inline Cmpt invariantI(const SymmTensor<Cmpt>& st)
{
    return tr(st);
}


//- Return the 2nd invariant of a SymmTensor
template<class Cmpt>
inline Cmpt invariantII(const SymmTensor<Cmpt>& st)
{
    return
    (
        st.xx()*st.yy() + st.yy()*st.zz() + st.xx()*st.zz()
      - st.xy()*st.xy() - st.yz()*st.yz() - st.xz()*st.xz()
    );
}


//- Return the 3rd invariant of a SymmTensor
template<class Cmpt>
inline Cmpt invariantIII(const SymmTensor<Cmpt>& st)
{
    return det(st);
}


//- Return the inner-product of a SymmTensor with itself
template<class Cmpt>
inline SymmTensor<Cmpt>
innerSqr(const SymmTensor<Cmpt>& st)
{
    return SymmTensor<Cmpt>
    (
        st.xx()*st.xx() + st.xy()*st.xy() + st.xz()*st.xz(),
        st.xx()*st.xy() + st.xy()*st.yy() + st.xz()*st.yz(),
        st.xx()*st.xz() + st.xy()*st.yz() + st.xz()*st.zz(),

        st.xy()*st.xy() + st.yy()*st.yy() + st.yz()*st.yz(),
        st.xy()*st.xz() + st.yy()*st.yz() + st.yz()*st.zz(),

        st.xz()*st.xz() + st.yz()*st.yz() + st.zz()*st.zz()
    );
}


//- Return the square of Frobenius norm of a SymmTensor as a Cmpt
template<class Cmpt>
inline Cmpt magSqr(const SymmTensor<Cmpt>& st)
{
    return Cmpt
    (
        mag(st.xx()*st.xx()) + 2*mag(st.xy()*st.xy()) + 2*mag(st.xz()*st.xz())
                             +   mag(st.yy()*st.yy()) + 2*mag(st.yz()*st.yz())
                                                      +   mag(st.zz()*st.zz())
    );
}


//- Return the square of a Vector as a SymmTensor
template<class Cmpt>
inline SymmTensor<Cmpt> sqr(const Vector<Cmpt>& v)
{
    return SymmTensor<Cmpt>
    (
        v.x()*v.x(), v.x()*v.y(), v.x()*v.z(),
                     v.y()*v.y(), v.y()*v.z(),
                                  v.z()*v.z()
    );
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Sum of a SphericalTensor and a SymmTensor
template<class Cmpt>
inline SymmTensor<Cmpt>
operator+(const SphericalTensor<Cmpt>& spt1, const SymmTensor<Cmpt>& st2)
{
    return SymmTensor<Cmpt>
    (
        spt1.ii() + st2.xx(), st2.xy(),             st2.xz(),
                              spt1.ii() + st2.yy(), st2.yz(),
                                                    spt1.ii() + st2.zz()
    );
}


//- Sum of a SymmTensor and a SphericalTensor
template<class Cmpt>
inline SymmTensor<Cmpt>
operator+(const SymmTensor<Cmpt>& st1, const SphericalTensor<Cmpt>& spt2)
{
    return SymmTensor<Cmpt>
    (
        st1.xx() + spt2.ii(), st1.xy(),             st1.xz(),
                              st1.yy() + spt2.ii(), st1.yz(),
                                                    st1.zz() + spt2.ii()
    );
}


//- Subtract a SymmTensor from a SphericalTensor
template<class Cmpt>
inline SymmTensor<Cmpt>
operator-(const SphericalTensor<Cmpt>& spt1, const SymmTensor<Cmpt>& st2)
{
    return SymmTensor<Cmpt>
    (
        spt1.ii() - st2.xx(), -st2.xy(),             -st2.xz(),
                               spt1.ii() - st2.yy(), -st2.yz(),
                                                      spt1.ii() - st2.zz()
    );
}


//- Subtract a SphericalTensor from a SymmTensor
template<class Cmpt>
inline SymmTensor<Cmpt>
operator-(const SymmTensor<Cmpt>& st1, const SphericalTensor<Cmpt>& spt2)
{
    return SymmTensor<Cmpt>
    (
        st1.xx() - spt2.ii(), st1.xy(),             st1.xz(),
                              st1.yy() - spt2.ii(), st1.yz(),
                                                    st1.zz() - spt2.ii()
    );
}


//- Return the Hodge dual of a SymmTensor as a Vector
template<class Cmpt>
inline Vector<Cmpt> operator*(const SymmTensor<Cmpt>& st)
{
    return Vector<Cmpt>(st.yz(), -st.xz(), st.xy());
}


//- Division of a SymmTensor by a Cmpt
template<class Cmpt>
inline SymmTensor<Cmpt>
operator/(const SymmTensor<Cmpt>& st, const Cmpt s)
{
    return SymmTensor<Cmpt>
    (
        st.xx()/s, st.xy()/s, st.xz()/s,
                   st.yy()/s, st.yz()/s,
                              st.zz()/s
    );
}


//- Inner-product of a SymmTensor and a SymmTensor
template<class Cmpt>
inline Tensor<Cmpt>
operator&(const SymmTensor<Cmpt>& st1, const SymmTensor<Cmpt>& st2)
{
    return Tensor<Cmpt>
    (
        st1.xx()*st2.xx() + st1.xy()*st2.xy() + st1.xz()*st2.xz(),
        st1.xx()*st2.xy() + st1.xy()*st2.yy() + st1.xz()*st2.yz(),
        st1.xx()*st2.xz() + st1.xy()*st2.yz() + st1.xz()*st2.zz(),

        st1.xy()*st2.xx() + st1.yy()*st2.xy() + st1.yz()*st2.xz(),
        st1.xy()*st2.xy() + st1.yy()*st2.yy() + st1.yz()*st2.yz(),
        st1.xy()*st2.xz() + st1.yy()*st2.yz() + st1.yz()*st2.zz(),

        st1.xz()*st2.xx() + st1.yz()*st2.xy() + st1.zz()*st2.xz(),
        st1.xz()*st2.xy() + st1.yz()*st2.yy() + st1.zz()*st2.yz(),
        st1.xz()*st2.xz() + st1.yz()*st2.yz() + st1.zz()*st2.zz()
    );
}


//- Inner-product of a SphericalTensor and a SymmTensor
template<class Cmpt>
inline SymmTensor<Cmpt>
operator&(const SphericalTensor<Cmpt>& spt1, const SymmTensor<Cmpt>& st2)
{
    return SymmTensor<Cmpt>
    (
        spt1.ii()*st2.xx(), spt1.ii()*st2.xy(), spt1.ii()*st2.xz(),
                            spt1.ii()*st2.yy(), spt1.ii()*st2.yz(),
                                                spt1.ii()*st2.zz()
    );
}


//- Inner-product of a SymmTensor and a SphericalTensor
template<class Cmpt>
inline SymmTensor<Cmpt>
operator&(const SymmTensor<Cmpt>& st1, const SphericalTensor<Cmpt>& spt2)
{
    return SymmTensor<Cmpt>
    (
        st1.xx()*spt2.ii(), st1.xy()*spt2.ii(), st1.xz()*spt2.ii(),
                            st1.yy()*spt2.ii(), st1.yz()*spt2.ii(),
                                                st1.zz()*spt2.ii()
    );
}


//- Inner-product of a SymmTensor and a Vector
template<class Cmpt>
inline Vector<Cmpt>
operator&(const SymmTensor<Cmpt>& st, const Vector<Cmpt>& v)
{
    return Vector<Cmpt>
    (
        st.xx()*v.x() + st.xy()*v.y() + st.xz()*v.z(),
        st.xy()*v.x() + st.yy()*v.y() + st.yz()*v.z(),
        st.xz()*v.x() + st.yz()*v.y() + st.zz()*v.z()
    );
}


//- Inner-product of a Vector and a SymmTensor
template<class Cmpt>
inline Vector<Cmpt>
operator&(const Vector<Cmpt>& v, const SymmTensor<Cmpt>& st)
{
    return Vector<Cmpt>
    (
        v.x()*st.xx() + v.y()*st.xy() + v.z()*st.xz(),
        v.x()*st.xy() + v.y()*st.yy() + v.z()*st.yz(),
        v.x()*st.xz() + v.y()*st.yz() + v.z()*st.zz()
    );
}


//- Double-inner-product of a SymmTensor and a SymmTensor
template<class Cmpt>
inline Cmpt
operator&&(const SymmTensor<Cmpt>& st1, const SymmTensor<Cmpt>& st2)
{
    return
    (
        st1.xx()*st2.xx() + 2*st1.xy()*st2.xy() + 2*st1.xz()*st2.xz()
                          +   st1.yy()*st2.yy() + 2*st1.yz()*st2.yz()
                                                +   st1.zz()*st2.zz()
    );
}


//- Double-inner-product of a SphericalTensor and a SymmTensor
template<class Cmpt>
inline Cmpt
operator&&(const SphericalTensor<Cmpt>& spt1, const SymmTensor<Cmpt>& st2)
{
    return (spt1.ii()*st2.xx() + spt1.ii()*st2.yy() + spt1.ii()*st2.zz());
}


//- Double-inner-product of a SymmTensor and a SphericalTensor
template<class Cmpt>
inline Cmpt
operator&&(const SymmTensor<Cmpt>& st1, const SphericalTensor<Cmpt>& spt2)
{
    return (st1.xx()*spt2.ii() + st1.yy()*spt2.ii() + st1.zz()*spt2.ii());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Cmpt>
class outerProduct<SymmTensor<Cmpt>, Cmpt>
{
public:

    typedef SymmTensor<Cmpt> type;
};

template<class Cmpt>
class outerProduct<Cmpt, SymmTensor<Cmpt>>
{
public:

    typedef SymmTensor<Cmpt> type;
};

template<class Cmpt>
class innerProduct<SymmTensor<Cmpt>, SymmTensor<Cmpt>>
{
public:

    typedef Tensor<Cmpt> type;
};

template<class Cmpt>
class innerProduct<SymmTensor<Cmpt>, Vector<Cmpt>>
{
public:

    typedef Vector<Cmpt> type;
};

template<class Cmpt>
class innerProduct<Vector<Cmpt>, SymmTensor<Cmpt>>
{
public:

    typedef Vector<Cmpt> type;
};


template<class Cmpt>
class typeOfSum<SphericalTensor<Cmpt>, SymmTensor<Cmpt>>
{
public:

    typedef SymmTensor<Cmpt> type;
};

template<class Cmpt>
class typeOfSum<SymmTensor<Cmpt>, SphericalTensor<Cmpt>>
{
public:

    typedef SymmTensor<Cmpt> type;
};

template<class Cmpt>
class innerProduct<SphericalTensor<Cmpt>, SymmTensor<Cmpt>>
{
public:

    typedef SymmTensor<Cmpt> type;
};

template<class Cmpt>
class innerProduct<SymmTensor<Cmpt>, SphericalTensor<Cmpt>>
{
public:

    typedef SymmTensor<Cmpt> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
