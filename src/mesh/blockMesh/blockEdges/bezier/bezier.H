/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenCFD Ltd.
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

Class
    Foam::blockEdges::bezier

Description
    Nth order bezier curve edge. Only interior control points should be
    specified. The outer points are taken as start and end. Note that the
    calculation of each point takes 0(N^2) time, where N is the number of
    control points. This edge type shouldn't therefore be used for finely
    discretised line data; polyLine or similar will be more appropriate for
    such cases. Beziers are useful for simple curved shapes such as aerofoils,
    or when you want an edge to match a specific direction at one, or both, or
    its endpoints. In comparison with BSplines, the grading of bezier edges
    should be smoother, and the code is much simpler. The algorithmic order is
    worse, however, and the edge will not follow the control points as closely.

SourceFiles
    bezier.C

\*---------------------------------------------------------------------------*/

#ifndef blockEdges_bezier_H
#define blockEdges_bezier_H

#include "blockEdge.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blockEdges
{

/*---------------------------------------------------------------------------*\
                           Class bezier Declaration
\*---------------------------------------------------------------------------*/

class bezier
:
    public blockEdge
{
private:

        //- Control points
        pointField control_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        bezier(const bezier&);

        //- Disallow default bitwise assignment
        void operator=(const bezier&);


public:

    // Static data members
    TypeName("bezier");


    // Constructors

        //- Construct from components
        bezier
        (
            const pointField& points,
            const label start,
            const label end,
            const pointField& control
        );

        //- Construct from Istream
        bezier
        (
            const dictionary& dict,
            const label index,
            const searchableSurfaces& geometry,
            const pointField&,
            Istream&
        );


        //- Destructor
        virtual ~bezier();


    // Member Functions

        //- Return the point position corresponding to the curve parameter
        //  0 <= lambda <= 1
        point position(const scalar lambda) const;

        //- Return the length of the curve
        scalar length() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End of namespace blockEdges
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
