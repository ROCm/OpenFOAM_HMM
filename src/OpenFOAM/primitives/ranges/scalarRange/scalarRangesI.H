/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

inline Foam::scalarRanges::scalarRanges(const std::string& str, bool report)
:
    List<scalarRange>(scalarRanges::parse(str, report))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::scalarRanges::contains(const scalar value) const
{
    return match(value);
}


inline bool Foam::scalarRanges::match(const scalar value) const
{
    for (const scalarRange& range : *this)
    {
        if (range.match(value))
        {
            return true;
        }
    }

    return false;
}


// ************************************************************************* //
