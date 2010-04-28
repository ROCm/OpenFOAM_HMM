/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::vectorTensorTransform::vectorTensorTransform()
:
    t_(vector::zero),
    R_(sphericalTensor::I)
{}

inline Foam::vectorTensorTransform::vectorTensorTransform
(
    const vector& t,
    const tensor& R
)
:
    t_(t),
    R_(R)
{}

inline Foam::vectorTensorTransform::vectorTensorTransform(const vector& t)
:
    t_(t),
    R_(sphericalTensor::I)
{}

inline Foam::vectorTensorTransform::vectorTensorTransform(const tensor& R)
:
    t_(vector::zero),
    R_(R)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::vector& Foam::vectorTensorTransform::t() const
{
    return t_;
}


inline const Foam::tensor& Foam::vectorTensorTransform::R() const
{
    return R_;
}


inline Foam::vector& Foam::vectorTensorTransform::t()
{
    return t_;
}


inline Foam::tensor& Foam::vectorTensorTransform::R()
{
    return R_;
}


inline Foam::vector Foam::vectorTensorTransform::transform
(
    const vector& v
) const
{
    return t() + (R() & v);
}


inline Foam::vector Foam::vectorTensorTransform::invTransform
(
    const vector& v
) const
{
    return (R().T() & (v - t()));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::vectorTensorTransform::operator=
(
    const vectorTensorTransform& tr
)
{
    t_ = tr.t_;
    R_ = tr.R_;
}


inline void Foam::vectorTensorTransform::operator&=
(
    const vectorTensorTransform& tr
)
{
    t_ += tr.t_;
    R_ = tr.R_ & R_;
}


inline void Foam::vectorTensorTransform::operator=(const vector& t)
{
    t_ = t;
}


inline void Foam::vectorTensorTransform::operator+=(const vector& t)
{
    t_ += t;
}


inline void Foam::vectorTensorTransform::operator-=(const vector& t)
{
    t_ -= t;
}


inline void Foam::vectorTensorTransform::operator=(const tensor& R)
{
    R_ = R;
}


inline void Foam::vectorTensorTransform::operator&=(const tensor& R)
{
    R_ = R & R_;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

inline bool Foam::operator==
(
    const vectorTensorTransform& tr1,
    const vectorTensorTransform& tr2
)
{
    return (tr1.t() == tr2.t() && tr1.R() == tr2.R());
}


inline bool Foam::operator!=
(
    const vectorTensorTransform& tr1,
    const vectorTensorTransform& tr2
)
{
    return !operator==(tr1, tr2);
}


inline Foam::vectorTensorTransform Foam::operator+
(
    const vectorTensorTransform& tr,
    const vector& t
)
{
    return vectorTensorTransform(tr.t() + t, tr.R());
}


inline Foam::vectorTensorTransform Foam::operator+
(
    const vector& t,
    const vectorTensorTransform& tr
)
{
    return vectorTensorTransform(t + tr.t(), tr.R());
}


inline Foam::vectorTensorTransform Foam::operator-
(
    const vectorTensorTransform& tr,
    const vector& t
)
{
    return vectorTensorTransform(tr.t() - t, tr.R());
}


inline Foam::vectorTensorTransform Foam::operator&
(
    const vectorTensorTransform& tr1,
    const vectorTensorTransform& tr2
)
{
    return vectorTensorTransform
    (
        tr1.t() + tr2.t(),
        tr1.R() & tr2.R()
    );
}


// ************************************************************************* //
