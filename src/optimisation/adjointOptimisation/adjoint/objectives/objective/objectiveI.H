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

inline const Foam::word& Foam::objective::objectiveName() const
{
    return objectiveName_;
}


inline bool Foam::objective::hasdJdb() const
{
    return bool(dJdbPtr_);
}


inline bool Foam::objective::hasBoundarydJdb() const
{
    return bool(bdJdbPtr_);
}


inline bool Foam::objective::hasdSdbMult() const
{
    return bool(bdSdbMultPtr_);
}


inline bool Foam::objective::hasdndbMult() const
{
    return bool(bdndbMultPtr_);
}


inline bool Foam::objective::hasdxdbMult() const
{
    return bool(bdxdbMultPtr_);
}


inline bool Foam::objective::hasdxdbDirectMult() const
{
    return bool(bdxdbDirectMultPtr_);
}


inline bool Foam::objective::hasBoundaryEdgeContribution() const
{
    return bool(bEdgeContribution_);
}


inline bool Foam::objective::hasDivDxDbMult() const
{
    return bool(divDxDbMultPtr_);
}


inline bool Foam::objective::hasGradDxDbMult() const
{
    return bool(gradDxDbMultPtr_);
}


inline bool Foam::objective::hasBoundarydJdStress() const
{
    return bool(bdJdStressPtr_);
}


inline bool Foam::objective::hasIntegrationStartTime() const
{
    return bool(integrationStartTimePtr_);
}


inline bool Foam::objective::hasIntegrationEndTime() const
{
    return bool(integrationEndTimePtr_);
}

// ************************************************************************* //
