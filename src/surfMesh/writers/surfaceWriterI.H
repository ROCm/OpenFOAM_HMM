/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

inline Foam::label Foam::surfaceWriter::nFields() const
{
    return nFields_;
}


inline Foam::label Foam::surfaceWriter::nFields(const label n)
{
    label old(nFields_);
    nFields_ = n;
    return old;
}


inline bool Foam::surfaceWriter::isPointData() const
{
    return isPointData_;
}


inline bool Foam::surfaceWriter::isPointData(const bool on)
{
    bool old(isPointData_);
    isPointData_ = on;
    return old;
}


inline bool Foam::surfaceWriter::useTimeDir() const
{
    return useTimeDir_;
}


inline bool Foam::surfaceWriter::useTimeDir(const bool on)
{
    bool old(useTimeDir_);
    useTimeDir_ = on;
    return old;
}


inline bool Foam::surfaceWriter::verbose() const
{
    return verbose_;
}


inline bool Foam::surfaceWriter::verbose(const bool on)
{
    bool old(verbose_);
    verbose_ = on;
    return old;
}


inline Foam::scalar Foam::surfaceWriter::mergeDim() const
{
    return mergeDim_;
}


inline Foam::scalar Foam::surfaceWriter::mergeDim(const scalar dist)
{
    scalar old(mergeDim_);
    mergeDim_ = dist;
    return old;
}


inline bool Foam::surfaceWriter::hasTime() const
{
    return currTime_.name().size();
}


inline const Foam::word& Foam::surfaceWriter::timeName() const
{
    return currTime_.name();
}


inline Foam::scalar Foam::surfaceWriter::timeValue() const
{
    return currTime_.name().empty() ? 0 : currTime_.value();
}


// ************************************************************************* //
