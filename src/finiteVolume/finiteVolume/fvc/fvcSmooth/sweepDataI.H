/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

#include "transform.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline bool Foam::sweepData::update
(
    const sweepData::sweepData& svf,
    const point& position,
    const scalar tol
)
{
    if (!valid())
    {
        operator=(svf);
        return true;
    }

    scalar myDist2 = magSqr(position - origin());

    if (myDist2 < SMALL)
    {
        if (svf.value() > value())
        {
            operator=(svf);
            return true;
        }
        else
        {
            return false;
        }
    }

    scalar dist2 = magSqr(position - svf.origin());

    if (dist2 < myDist2)
    {
        operator=(svf);
        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::sweepData::sweepData()
:
    value_(-GREAT),
    origin_(vector::max)
{}


inline Foam::sweepData::sweepData(const scalar value, const point& origin)
:
    value_(value),
    origin_(origin)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::sweepData::valid() const
{
    return value_ > -SMALL;
}


inline bool Foam::sweepData::sameGeometry
(
    const polyMesh&,
    const sweepData&,
    const scalar
) const
{
    return true;
}


inline void Foam::sweepData::leaveDomain
(
    const polyMesh&,
    const polyPatch&,
    const label,
    const point& faceCentre
)
{
    origin_ -= faceCentre;
}


inline void Foam::sweepData::transform
(
    const polyMesh&,
    const tensor& rotTensor
)
{
    origin_ = Foam::transform(rotTensor, origin_);
}


inline void Foam::sweepData::enterDomain
(
    const polyMesh&,
    const polyPatch&,
    const label,
    const point& faceCentre
)
{
    // back to absolute form
    origin_ += faceCentre;
}


inline bool Foam::sweepData::updateCell
(
    const polyMesh& mesh,
    const label thisCellI,
    const label,
    const sweepData& svf,
    const scalar tol
)
{
    return update(svf, mesh.cellCentres()[thisCellI], tol);
}


inline bool Foam::sweepData::updateFace
(
    const polyMesh& mesh,
    const label thisFaceI,
    const label,
    const sweepData& svf,
    const scalar tol
)
{
    return update(svf, mesh.faceCentres()[thisFaceI], tol);
}


// Update this (face) with coupled face information.
inline bool Foam::sweepData::updateFace
(
    const polyMesh& mesh,
    const label thisFaceI,
    const sweepData& svf,
    const scalar tol
)
{
    return update(svf, mesh.faceCentres()[thisFaceI], tol);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::sweepData::operator=
(
    const scalar value
)
{
    value_ = value;
}


inline bool Foam::sweepData::operator==
(
    const sweepData& rhs
) const
{
    return origin() == rhs.origin();
}


inline bool Foam::sweepData::operator!=
(
    const sweepData& rhs
) const
{
    return !(*this == rhs);
}


// ************************************************************************* //
