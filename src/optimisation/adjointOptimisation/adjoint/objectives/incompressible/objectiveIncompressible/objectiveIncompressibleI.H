/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
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

inline bool Foam::objectiveIncompressible::hasdJdv() const
{
    return bool(dJdvPtr_);
}


inline bool Foam::objectiveIncompressible::hasdJdp() const
{
    return bool(dJdpPtr_);
}


inline bool Foam::objectiveIncompressible::hasdJdT() const
{
    return bool(dJdTPtr_);
}


inline bool Foam::objectiveIncompressible::hasdJdTMVar1() const
{
    return bool(dJdTMvar1Ptr_);
}


inline bool Foam::objectiveIncompressible::hasdJdTMVar2() const
{
    return bool(dJdTMvar2Ptr_);
}


inline bool Foam::objectiveIncompressible::hasBoundarydJdv() const
{
    return bool(bdJdvPtr_);
}


inline bool Foam::objectiveIncompressible::hasBoundarydJdvn() const
{
    return bool(bdJdvnPtr_);
}


inline bool Foam::objectiveIncompressible::hasBoundarydJdvt() const
{
    return bool(bdJdvtPtr_);
}


inline bool Foam::objectiveIncompressible::hasBoundarydJdp() const
{
    return bool(bdJdpPtr_);
}


inline bool Foam::objectiveIncompressible::hasBoundarydJdT() const
{
    return bool(bdJdTPtr_);
}


inline bool Foam::objectiveIncompressible::hasBoundarydJdTMVar1() const
{
    return bool(bdJdTMvar1Ptr_);
}


inline bool Foam::objectiveIncompressible::hasBoundarydJdTMVar2() const
{
    return bool(bdJdTMvar2Ptr_);
}


// ************************************************************************* //
