/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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


inline bool Foam::lumpedPointMovement::empty() const
{
    return locations_.empty();
}


inline Foam::label Foam::lumpedPointMovement::size() const
{
    return locations_.size();
}


inline const Foam::vector& Foam::lumpedPointMovement::axis() const
{
    return axis_;
}


inline const Foam::scalarField& Foam::lumpedPointMovement::locations() const
{
    return locations_;
}


inline Foam::scalar Foam::lumpedPointMovement::division() const
{
    return division_;
}


inline Foam::label Foam::lumpedPointMovement::ownerId() const
{
    return ownerId_;
}


inline void Foam::lumpedPointMovement::ownerId(label id)
{
    ownerId_ = id;
}


inline const Foam::scalarField& Foam::lumpedPointMovement::thresholds() const
{
    if (!thresholdPtr_)
    {
        calcThresholds();
    }

    return *thresholdPtr_;
}


inline Foam::label
Foam::lumpedPointMovement::threshold(const point& position) const
{
    return threshold(position & axis_);
}


inline const Foam::externalFileCoupler&
Foam::lumpedPointMovement::coupler() const
{
    return coupler_;
}


inline Foam::externalFileCoupler& Foam::lumpedPointMovement::coupler()
{
    return coupler_;
}


//- The initial state (positions/rotations)
inline const Foam::lumpedPointState& Foam::lumpedPointMovement::state0() const
{
    return state0_;
}


inline const Foam::lumpedPointState& Foam::lumpedPointMovement::state() const
{
    return state_;
}


inline Foam::scalar Foam::lumpedPointMovement::relax() const
{
    return relax_;
}


inline Foam::scalar& Foam::lumpedPointMovement::relax()
{
    return relax_;
}


inline const Foam::word& Foam::lumpedPointMovement::inputName() const
{
    return inputName_;
}


inline const Foam::word& Foam::lumpedPointMovement::outputName() const
{
    return outputName_;
}


inline Foam::lumpedPointState::inputFormatType
Foam::lumpedPointMovement::inputFormat() const
{
    return inputFormat_;
}


Foam::lumpedPointMovement::outputFormatType
Foam::lumpedPointMovement::outputFormat() const
{
    return outputFormat_;
}


inline bool Foam::lumpedPointMovement::hasMapping() const
{
    return !faceZones_.empty();
}


inline const Foam::List<Foam::labelList>&
Foam::lumpedPointMovement::zones() const
{
    return faceZones_;
}


// ************************************************************************* //
