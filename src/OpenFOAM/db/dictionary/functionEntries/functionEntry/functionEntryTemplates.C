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

#include "functionEntry.H"
#include "ISstream.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class StringType>
Foam::List<StringType>
Foam::functionEntry::readStringList(Istream& is)
{
    ISstream& iss = dynamic_cast<ISstream&>(is);
    token firstToken(iss);

    List<StringType> list;

    if (firstToken.isWord() || firstToken.isString())
    {
        // The first token appears viable as non-list
        // - treated like list with one entry

        iss.putBack(firstToken);

        list.setSize(1);

        iss >> list[0];
    }
    else
    {
        iss.putBack(firstToken);
        iss >> list;
    }

    return list;
}


// ************************************************************************* //
