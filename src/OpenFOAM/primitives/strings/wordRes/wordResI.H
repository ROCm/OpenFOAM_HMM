/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

inline Foam::wordRes::wordRes
(
    const UList<wordRe>& list
)
:
    list_(list)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::wordRes::size() const
{
    return list_.size();
}


inline bool Foam::wordRes::empty() const
{
    return list_.empty();
}


inline bool Foam::wordRes::match
(
    const std::string& text,
    const bool literal
) const
{
    const label n = list_.size();

    for (label i = 0; i < n; ++i)
    {
        if (list_[i].match(text, literal))
        {
            return true;
        }
    }

    return false;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline const Foam::UList<Foam::wordRe>& Foam::wordRes::operator()() const
{
    return list_;
}


inline bool Foam::wordRes::operator()(const std::string& text) const
{
    return match(text);
}


inline const Foam::wordRe& Foam::wordRes::operator[](const label i) const
{
    return list_[i];
}


// ************************************************************************* //
