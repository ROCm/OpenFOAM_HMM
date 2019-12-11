/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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
    return dJdvPtr_.valid();
}


inline bool Foam::objectiveIncompressible::hasdJdp() const
{
    return dJdpPtr_.valid();
}


inline bool Foam::objectiveIncompressible::hasdJdT() const
{
    return dJdTPtr_.valid();
}


inline bool Foam::objectiveIncompressible::hasdJdTMVar1() const
{
    return dJdTMvar1Ptr_.valid();
}


inline bool Foam::objectiveIncompressible::hasdJdTMVar2() const
{
    return dJdTMvar2Ptr_.valid();
}


inline bool Foam::objectiveIncompressible::hasBoundarydJdv() const
{
    return bdJdvPtr_.valid();
}


inline bool Foam::objectiveIncompressible::hasBoundarydJdvn() const
{
    return bdJdvnPtr_.valid();
}


inline bool Foam::objectiveIncompressible::hasBoundarydJdvt() const
{
    return bdJdvtPtr_.valid();
}


inline bool Foam::objectiveIncompressible::hasBoundarydJdp() const
{
    return bdJdpPtr_.valid();
}


inline bool Foam::objectiveIncompressible::hasBoundarydJdT() const
{
    return bdJdTPtr_.valid();
}


inline bool Foam::objectiveIncompressible::hasBoundarydJdTMVar1() const
{
    return bdJdTMvar1Ptr_.valid();
}


inline bool Foam::objectiveIncompressible::hasBoundarydJdTMVar2() const
{
    return bdJdTMvar2Ptr_.valid();
}


// ************************************************************************* //
