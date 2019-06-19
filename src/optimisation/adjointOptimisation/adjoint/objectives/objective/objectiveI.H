/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2007-2019 PCOpt/NTUA
                            | Copyright (C) 2013-2019 FOSS GP
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


inline bool Foam::objective::hasdJdb()
{
    return dJdbPtr_.valid();
}


inline bool Foam::objective::hasBoundarydJdb()
{
    return bdJdbPtr_.valid();
}


inline bool Foam::objective::hasdSdbMult()
{
    return bdSdbMultPtr_.valid();
}


inline bool Foam::objective::hasdndbMult()
{
    return bdndbMultPtr_.valid();
}


inline bool Foam::objective::hasdxdbMult()
{
    return bdxdbMultPtr_.valid();
}


inline bool Foam::objective::hasdxdbDirectMult()
{
    return bdxdbDirectMultPtr_.valid();
}


inline bool Foam::objective::hasBoundaryEdgeContribution()
{
    return bEdgeContribution_.valid();
}


inline bool Foam::objective::hasDivDxDbMult()
{
    return divDxDbMultPtr_.valid();
}


inline bool Foam::objective::hasGradDxDbMult()
{
    return gradDxDbMultPtr_.valid();
}


inline bool Foam::objective::hasBoundarydJdStress()
{
    return bdJdStressPtr_.valid();
}


// ************************************************************************* //
