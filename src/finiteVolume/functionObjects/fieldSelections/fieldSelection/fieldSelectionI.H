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

inline Foam::HashSet<Foam::wordRe>
Foam::functionObjects::fieldSelection::filters() const
{
    HashSet<wordRe> f;
    for (const fieldInfo& fi : *this)
    {
        f.insert(fi.name());
    }

    return f;
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

inline const Foam::List<Foam::functionObjects::fieldInfo>&
Foam::functionObjects::fieldSelection::selection() const
{
    return selection_;
}


inline Foam::wordHashSet
Foam::functionObjects::fieldSelection::selectionNames() const
{
    wordHashSet names;
    for (const fieldInfo& fi : selection_)
    {
        names.insert(fi.name());
    }

    return names;
}


// ************************************************************************* //
