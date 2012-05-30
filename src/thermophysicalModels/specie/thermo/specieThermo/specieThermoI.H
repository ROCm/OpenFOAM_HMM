/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "specieThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline Foam::specieThermo<Thermo, Type>::specieThermo
(
    const Thermo& sp
)
:
    Thermo(sp)
{}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::specieThermo<Thermo, Type>::T
(
    scalar f,
    scalar T0,
    scalar (specieThermo<Thermo, Type>::*F)(const scalar) const,
    scalar (specieThermo<Thermo, Type>::*dFdT)(const scalar) const,
    scalar (specieThermo<Thermo, Type>::*limit)(const scalar) const
) const
{
    scalar Test = T0;
    scalar Tnew = T0;
    scalar Ttol = T0*tol_;
    int    iter = 0;

    do
    {
        Test = Tnew;
        Tnew =
            (this->*limit)(Test - ((this->*F)(Test) - f)/(this->*dFdT)(Test));

        if (iter++ > maxIter_)
        {
            FatalErrorIn
            (
                "specieThermo<Thermo, Type>::T(scalar f, scalar T0, "
                "scalar (specieThermo<Thermo, Type>::*F)"
                "(const scalar) const, "
                "scalar (specieThermo<Thermo, Type>::*dFdT)"
                "(const scalar) const"
                ") const"
            )   << "Maximum number of iterations exceeded"
                << abort(FatalError);
        }

    } while (mag(Tnew - Test) > Ttol);

    return Tnew;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline Foam::specieThermo<Thermo, Type>::specieThermo
(
    const word& name,
    const specieThermo& st
)
:
    Thermo(name, st)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline Foam::word
Foam::specieThermo<Thermo, Type>::heName()
{
    return Type<specieThermo<Thermo, Type> >::name();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::specieThermo<Thermo, Type>::he(const scalar T) const
{
    return Type<specieThermo<Thermo, Type> >::he(*this, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::specieThermo<Thermo, Type>::cv(const scalar T) const
{
    return this->cp(T) - this->RR;
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::specieThermo<Thermo, Type>::cpv(const scalar T) const
{
    return Type<specieThermo<Thermo, Type> >::cpv(*this, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::specieThermo<Thermo, Type>::gamma(const scalar T) const
{
    scalar CP = this->cp(T);
    return CP/(CP - this->RR);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::specieThermo<Thermo, Type>::cpBycpv(const scalar T) const
{
    return Type<specieThermo<Thermo, Type> >::cpBycpv(*this, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::specieThermo<Thermo, Type>::es(const scalar T) const
{
    return this->hs(T) - this->RR*(T - this->Tstd);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::specieThermo<Thermo, Type>::ea(const scalar T) const
{
    return this->ha(T) - this->RR*(T - this->Tstd);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::specieThermo<Thermo, Type>::g(const scalar T) const
{
    return this->ha(T) - T*this->s(T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::specieThermo<Thermo, Type>::a(const scalar T) const
{
    return this->ea(T) - T*this->s(T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::specieThermo<Thermo, Type>::Cpv(const scalar T) const
{
    return this->cpv(T)/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::specieThermo<Thermo, Type>::Cp(const scalar T) const
{
    return this->cp(T)/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::specieThermo<Thermo, Type>::Cv(const scalar T) const
{
    return this->cv(T)/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::specieThermo<Thermo, Type>::HE(const scalar T) const
{
    return Type<specieThermo<Thermo, Type> >::HE(*this, T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::specieThermo<Thermo, Type>::H(const scalar T) const
{
    return this->h(T)/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::specieThermo<Thermo, Type>::Hs(const scalar T) const
{
    return this->hs(T)/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::specieThermo<Thermo, Type>::Hc() const
{
    return this->hc()/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::specieThermo<Thermo, Type>::Ha(const scalar T) const
{
    return this->ha(T)/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::specieThermo<Thermo, Type>::S(const scalar T) const
{
    return this->s(T)/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::specieThermo<Thermo, Type>::E(const scalar T) const
{
    return this->e(T)/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::specieThermo<Thermo, Type>::Es(const scalar T) const
{
    return this->es(T)/this->W();
}

template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::specieThermo<Thermo, Type>::Ea(const scalar T) const
{
    return this->ea(T)/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::specieThermo<Thermo, Type>::G(const scalar T) const
{
    return this->g(T)/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::specieThermo<Thermo, Type>::A(const scalar T) const
{
    return this->a(T)/this->W();
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::specieThermo<Thermo, Type>::K(const scalar T) const
{
    scalar arg = -this->nMoles()*this->g(T)/(this->RR*T);

    if (arg < 600.0)
    {
        return ::exp(arg);
    }
    else
    {
        return VGREAT;
    }
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::specieThermo<Thermo, Type>::Kp(const scalar T) const
{
    return K(T);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar
Foam::specieThermo<Thermo, Type>::Kc(const scalar T) const
{
    if (equal(this->nMoles(), SMALL))
    {
        return Kp(T);
    }
    else
    {
        return Kp(T)*pow(this->Pstd/(this->RR*T), this->nMoles());
    }
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::specieThermo<Thermo, Type>::Kx
(
    const scalar T,
    const scalar p
) const
{
    if (equal(this->nMoles(), SMALL))
    {
        return Kp(T);
    }
    else
    {
        return Kp(T)*pow(this->Pstd/p, this->nMoles());
    }
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::specieThermo<Thermo, Type>::Kn
(
    const scalar T,
    const scalar p,
    const scalar n
) const
{
    if (equal(this->nMoles(), SMALL))
    {
        return Kp(T);
    }
    else
    {
        return Kp(T)*pow(n*this->Pstd/p, this->nMoles());
    }
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::specieThermo<Thermo, Type>::THE
(
    const scalar he,
    const scalar T0
) const
{
    return Type<specieThermo<Thermo, Type> >::THE(*this, he, T0);
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::specieThermo<Thermo, Type>::THs
(
    const scalar hs,
    const scalar T0
) const
{
    return T
    (
        hs,
        T0,
        &specieThermo<Thermo, Type>::Hs,
        &specieThermo<Thermo, Type>::Cp,
        &specieThermo<Thermo, Type>::limit
    );
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::specieThermo<Thermo, Type>::THa
(
    const scalar ht,
    const scalar T0
) const
{
    return T
    (
        ht,
        T0,
        &specieThermo<Thermo, Type>::Ha,
        &specieThermo<Thermo, Type>::Cp,
        &specieThermo<Thermo, Type>::limit
    );
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::specieThermo<Thermo, Type>::TEs
(
    const scalar es,
    const scalar T0
) const
{
    return T
    (
        es,
        T0,
        &specieThermo<Thermo, Type>::Es,
        &specieThermo<Thermo, Type>::Cv,
        &specieThermo<Thermo, Type>::limit
    );
}


template<class Thermo, template<class> class Type>
inline Foam::scalar Foam::specieThermo<Thermo, Type>::TEa
(
    const scalar e,
    const scalar T0
) const
{
    return T
    (
        ea,
        T0,
        &specieThermo<Thermo, Type>::Ea,
        &specieThermo<Thermo, Type>::Cv,
        &specieThermo<Thermo, Type>::limit
    );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline void Foam::specieThermo<Thermo, Type>::operator+=
(
    const specieThermo<Thermo, Type>& st
)
{
    Thermo::operator+=(st);
}


template<class Thermo, template<class> class Type>
inline void Foam::specieThermo<Thermo, Type>::operator-=
(
    const specieThermo<Thermo, Type>& st
)
{
    Thermo::operator-=(st);
}


template<class Thermo, template<class> class Type>
inline void Foam::specieThermo<Thermo, Type>::operator*=(const scalar s)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo, template<class> class Type>
inline Foam::specieThermo<Thermo, Type> Foam::operator+
(
    const specieThermo<Thermo, Type>& st1,
    const specieThermo<Thermo, Type>& st2
)
{
    return specieThermo<Thermo, Type>
    (
        static_cast<const Thermo&>(st1) + static_cast<const Thermo&>(st2)
    );
}


template<class Thermo, template<class> class Type>
inline Foam::specieThermo<Thermo, Type> Foam::operator-
(
    const specieThermo<Thermo, Type>& st1,
    const specieThermo<Thermo, Type>& st2
)
{
    return specieThermo<Thermo, Type>
    (
        static_cast<const Thermo&>(st1) - static_cast<const Thermo&>(st2)
    );
}


template<class Thermo, template<class> class Type>
inline Foam::specieThermo<Thermo, Type> Foam::operator*
(
    const scalar s,
    const specieThermo<Thermo, Type>& st
)
{
    return specieThermo<Thermo, Type>
    (
        s*static_cast<const Thermo&>(st)
    );
}


template<class Thermo, template<class> class Type>
inline Foam::specieThermo<Thermo, Type> Foam::operator==
(
    const specieThermo<Thermo, Type>& st1,
    const specieThermo<Thermo, Type>& st2
)
{
    return st2 - st1;
}


// ************************************************************************* //
