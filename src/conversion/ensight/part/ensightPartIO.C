/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

Description
    Output for ensightPart

\*---------------------------------------------------------------------------*/

#include "ensightPart.H"
#include "dictionary.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ensightPart::reconstruct(Istream& is)
{
    dictionary dict(is);
    dict.lookup("id") >> number_;
    dict.lookup("name") >> name_;

    offset_ = 0;
    dict.readIfPresent("offset", offset_);

    // populate elemLists_
    elemLists_.setSize(elementTypes().size());

    forAll(elementTypes(), elemI)
    {
        word key(elementTypes()[elemI]);

        elemLists_[elemI].clear();
        dict.readIfPresent(key, elemLists_[elemI]);

        size_ += elemLists_[elemI].size();
    }

    is.check("ensightPart::reconstruct(Istream&)");
}


bool Foam::ensightPart::writeSummary(Ostream& os) const
{
    os  << indent << type() << nl
        << indent << token::BEGIN_BLOCK << incrIndent << nl;

    // Ensight starts with 1
    os.writeKeyword("id") << (number() + 1) << token::END_STATEMENT << nl;
    os.writeKeyword("name") << name() << token::END_STATEMENT << nl;
    os.writeKeyword("offset") << offset() << token::END_STATEMENT << nl;
    os.writeKeyword("size") << size() << token::END_STATEMENT << nl;

    os  << decrIndent << indent << token::END_BLOCK << nl << endl;

    return true;
}


bool Foam::ensightPart::writeData(Ostream& os) const
{
    os  << indent << type() << nl
        << indent << token::BEGIN_BLOCK << incrIndent << nl;

    os.writeKeyword("id") << number() << token::END_STATEMENT << nl;
    os.writeKeyword("name") << name() << token::END_STATEMENT << nl;
    os.writeKeyword("offset") << offset() << token::END_STATEMENT << nl;

    forAll(elementTypes(), typeI)
    {
        word key(elementTypes()[typeI]);
        if (elemLists_[typeI].size())
        {
            elemLists_[typeI].writeEntry(key, os);
        }
    }

    os  << decrIndent << indent << token::END_BLOCK << nl << endl;

    return true;
}


void Foam::ensightPart::writeGeometry
(
    ensightGeoFile& os,
    const pointField& points
) const
{
    if (size())
    {
        const localPoints ptList = calcLocalPoints();
        const labelUList& pointMap = ptList.list;

        os.beginPart(number(), name());
        os.beginCoordinates(ptList.nPoints);

        for (direction cmpt=0; cmpt < point::nComponents; ++cmpt)
        {
            forAll(pointMap, ptI)
            {
                if (pointMap[ptI] > -1)
                {
                    os.write(points[ptI].component(cmpt));
                    os.newline();
                }
            }
        }

        // write parts
        forAll(elementTypes(), elemI)
        {
            if (elemLists_[elemI].size())
            {
                writeConnectivity
                (
                    os,
                    elementTypes()[elemI],
                    elemLists_[elemI],
                    pointMap
                );
            }
        }
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ensightPart& part
)
{
    part.writeData(os);
    return os;
}


Foam::ensightGeoFile& Foam::operator<<
(
    ensightGeoFile& os,
    const ensightPart& part
)
{
    part.writeGeometry(os);
    return os;
}


// ************************************************************************* //
