/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::fa::faceSetOption::timeStart() const
{
    return timeStart_;
}


inline Foam::scalar Foam::fa::faceSetOption::duration() const
{
    return duration_;
}


inline bool Foam::fa::faceSetOption::inTimeLimits(const scalar time) const
{
    return
    (
        (timeStart_ < 0)
     ||
        (
            (mesh_.time().value() >= timeStart_)
         && (mesh_.time().value() <= (timeStart_ + duration_))
        )
    );
}


inline const Foam::fa::faceSetOption::selectionModeType&
Foam::fa::faceSetOption::selectionMode() const
{
    return selectionMode_;
}


inline const Foam::word& Foam::fa::faceSetOption::faceSetName() const
{
    return faceSetName_;
}


inline Foam::scalar Foam::fa::faceSetOption::A() const
{
    return A_;
}


inline const Foam::labelList& Foam::fa::faceSetOption::faces() const
{
    return faces_;
}


inline Foam::scalar& Foam::fa::faceSetOption::timeStart()
{
    return timeStart_;
}


inline Foam::scalar& Foam::fa::faceSetOption::duration()
{
    return duration_;
}


// ************************************************************************* //
