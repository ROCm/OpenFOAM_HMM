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
    Foam::simpleSplineEdge

Description
    The actual access class for Bspline

SourceFiles
    simpleSplineEdge.C

\*---------------------------------------------------------------------------*/

#ifndef simpleSplineEdge_H
#define simpleSplineEdge_H

#include "curvedEdge.H"
#include "BSpline.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class simpleSplineEdge Declaration
\*---------------------------------------------------------------------------*/

class simpleSplineEdge
:
    public curvedEdge,
    public BSpline
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        simpleSplineEdge(const simpleSplineEdge&);

        //- Disallow default bitwise assignment
        void operator=(const simpleSplineEdge&);


public:

    //- Runtime type information
    TypeName("simpleSpline");


    // Constructors

        //- Construct from components
        simpleSplineEdge
        (
            const pointField&,
            const label start,
            const label end,
            const pointField& otherKnots
        );

        //- Construct from components
        simpleSplineEdge
        (
            const pointField&,
            const label start,
            const label end,
            const pointField& otherKnots,
            const vector& fstend,
            const vector& sndend
        );

        //- Construct from Istream setting pointsList
        simpleSplineEdge(const pointField&, Istream&);


    // Destructor

        virtual ~simpleSplineEdge()
        {}


    // Member Functions

        //- Return the point position corresponding to the curve parameter
        //  0 <= lambda <= 1
        point position(const scalar) const;

        //- Return the length of the simple spline curve
        scalar length() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
