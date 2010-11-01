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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline bool Foam::smoothData::update
(
    const smoothData& svf,
    const scalar scale,
    const scalar tol
)
{
    if (!valid() || (value_ < VSMALL))
    {
        // My value not set - take over neighbour
        value_ = svf.value()/scale;

        // Something changed - let caller know
        return true;
    }
    else if (svf.value() > (1 + tol)*scale*value_)
    {
        // Neighbour is too big for me - Up my value
        value_ = svf.value()/scale;

        // Something changed - let caller know
        return true;
    }
    else
    {
        // Neighbour is not too big for me or change is too small
        // Nothing changed
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::smoothData::smoothData()
:
    value_(-GREAT)
{}


inline Foam::smoothData::smoothData(const scalar value)
:
    value_(value)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::smoothData::valid() const
{
    return value_ > -SMALL;
}


inline bool Foam::smoothData::sameGeometry
(
    const polyMesh&,
    const smoothData&,
    const scalar
) const
{
    return true;
}


inline void Foam::smoothData::leaveDomain
(
    const polyMesh&,
    const polyPatch&,
    const label,
    const point&
)
{}


inline void Foam::smoothData::transform
(
    const polyMesh&,
    const tensor&
)
{}


inline void Foam::smoothData::enterDomain
(
    const polyMesh&,
    const polyPatch&,
    const label,
    const point&
)
{}


inline bool Foam::smoothData::updateCell
(
    const polyMesh&,
    const label,
    const label,
    const smoothData& svf,
    const scalar tol
)
{
    // Take over info from face if more than deltaRatio larger
    return update(svf, maxRatio, tol);
}


inline bool Foam::smoothData::updateFace
(
    const polyMesh&,
    const label,
    const label,
    const smoothData& svf,
    const scalar tol
)
{
    // Take over information from cell without any scaling (scale = 1.0)
    return update(svf, 1.0, tol);
}


// Update this (face) with coupled face information.
inline bool Foam::smoothData::updateFace
(
    const polyMesh&,
    const label,
    const smoothData& svf,
    const scalar tol
)
{
    // Take over information from coupled face without any scaling (scale = 1.0)
    return update(svf, 1.0, tol);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::smoothData::operator=
(
    const scalar value
)
{
    value_ = value;
}


inline bool Foam::smoothData::operator==
(
    const smoothData& rhs
) const
{
    return value_ == rhs.value();
}


inline bool Foam::smoothData::operator!=
(
    const smoothData& rhs
) const
{
    return !(*this == rhs);
}


// ************************************************************************* //
