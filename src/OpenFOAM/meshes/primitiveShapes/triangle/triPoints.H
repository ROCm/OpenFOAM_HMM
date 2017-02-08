/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

Class
    Foam::triPoints

Description
    Triangle storage. Null constructable (unfortunately triangle<point, point>
    is not)

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef triPoints_H
#define triPoints_H

#include "FixedList.H"
#include "treeBoundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class triPoints Declaration
\*---------------------------------------------------------------------------*/

class triPoints
:
    public FixedList<point, 3>
{
public:

    // Constructors

        //- Construct null
        inline triPoints()
        {}

        //- Construct from points
        inline triPoints
        (
            const point& a,
            const point& b,
            const point& c
        )
        {
            operator[](0) = a;
            operator[](1) = b;
            operator[](2) = c;
        }


    // Member Functions

        //- Calculate the bounding box
        inline treeBoundBox bounds() const
        {
            treeBoundBox bb(operator[](0));
            for (label i = 1; i < size(); ++i)
            {
                bb.add(operator[](i));
            }
            return bb;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
