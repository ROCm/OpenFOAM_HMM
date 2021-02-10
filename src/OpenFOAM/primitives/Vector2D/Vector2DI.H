/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Cmpt>
inline Foam::Vector2D<Cmpt>::Vector2D(const Foam::zero)
:
    Vector2D::vsType(Zero)
{}


template<class Cmpt>
inline Foam::Vector2D<Cmpt>::Vector2D
(
    const VectorSpace<Vector2D<Cmpt>, Cmpt, 2>& vs
)
:
    Vector2D::vsType(vs)
{}


template<class Cmpt>
inline Foam::Vector2D<Cmpt>::Vector2D(const Cmpt& vx, const Cmpt& vy)
{
    this->v_[X] = vx;
    this->v_[Y] = vy;
}


template<class Cmpt>
inline Foam::Vector2D<Cmpt>::Vector2D(Istream& is)
:
    Vector2D::vsType(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cmpt>
inline const Cmpt& Foam::Vector2D<Cmpt>::x() const
{
    return this->v_[X];
}


template<class Cmpt>
inline const Cmpt& Foam::Vector2D<Cmpt>::y() const
{
    return this->v_[Y];
}


template<class Cmpt>
inline Cmpt& Foam::Vector2D<Cmpt>::x()
{
    return this->v_[X];
}


template<class Cmpt>
inline Cmpt& Foam::Vector2D<Cmpt>::y()
{
    return this->v_[Y];
}


template<class Cmpt>
inline Foam::Vector2D<Cmpt>& Foam::Vector2D<Cmpt>::normalise()
{
    const scalar s(Foam::mag(*this));

    if (s < ROOTVSMALL)
    {
        *this = Zero;
    }
    else
    {
        *this /= s;
    }

    return *this;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Cmpt>
inline typename innerProduct<Vector2D<Cmpt>, Vector2D<Cmpt>>::type
operator&(const Vector2D<Cmpt>& v1, const Vector2D<Cmpt>& v2)
{
    return Cmpt(v1.x()*v2.x() + v1.y()*v2.y());
}


template<class Cmpt>
inline scalar Vector2D<Cmpt>::perp(const Vector2D<Cmpt>& b) const
{
    return x()*b.y() - y()*b.x();
}


template<class Cmpt>
inline bool Vector2D<Cmpt>::isClose
(
    const Vector2D<Cmpt>& b,
    const scalar tol
) const
{
    return (mag(x() - b.x()) < tol && mag(y() - b.y()) < tol);
}


template<class Cmpt>
inline Vector2D<Cmpt> operator*(const Cmpt& s, const Vector2D<Cmpt>& v)
{
    return Vector2D<Cmpt>(s*v.x(), s*v.y());
}


template<class Cmpt>
inline Vector2D<Cmpt> operator*(const Vector2D<Cmpt>& v, const Cmpt& s)
{
    return s*v;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
