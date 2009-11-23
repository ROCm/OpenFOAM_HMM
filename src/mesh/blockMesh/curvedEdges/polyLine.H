/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::polyLine

Description
    Define a series of control points, which can also be interpreted as a
    series of straight line segments.

    This is the basic polyLine class which implements just the line
    (no topology - it is not derived from curvedEdge)

SourceFiles
    polyLine.C

\*---------------------------------------------------------------------------*/

#ifndef polyLine_H
#define polyLine_H

#include "pointField.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class polyLine Declaration
\*---------------------------------------------------------------------------*/


class polyLine
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        polyLine(const polyLine&);

        //- Disallow default bitwise assignment
        void operator=(const polyLine&);

protected:

    // Protected data

        //- The control points or ends of each segmen
        pointField controlPoints_;

        //- The rational (0-1) cumulative distance for each control-point
        scalarList distances_;

        //- The real line length
        scalar lineLength_;

    // Protected member functions

        //- Precalculate the rational cumulative distances and the line-length
        void calcDistances();

public:

    // Constructors

        //- Construct from components
        polyLine(const pointField&);


    // Member Functions

        //- Return const-access to the control-points
        const pointField& controlPoints() const;

        //- Return the point position corresponding to the curve parameter
        //  0 <= lambda <= 1
        point position(const scalar) const;

        //- Return the length of the curve
        scalar length() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
