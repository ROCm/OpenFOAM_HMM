/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Class
    Foam::tetCell

Description
    A tetrahedral cell primitive.

    It is important that the ordering of edges is the same for a tetrahedron
    class, a tetrahedron cell shape model and a tetCell

SourceFiles
    tetCell.C
    tetCellI.H

\*---------------------------------------------------------------------------*/

#ifndef tetCell_H
#define tetCell_H

#include "FixedList.H"
#include "triFace.H"
#include "edge.H"
#include "pointField.H"
#include "tetPointRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class cellShape;

/*---------------------------------------------------------------------------*\
                           Class tetCell Declaration
\*---------------------------------------------------------------------------*/

class tetCell
:
    public FixedList<label, 4>
{
public:

    // Constructors

        //- Default construct, with invalid point labels (-1)
        inline tetCell();

        //- Construct from four point labels
        inline tetCell
        (
            const label a,
            const label b,
            const label c,
            const label d
        );

        //- Construct from an initializer list of four point labels
        inline explicit tetCell(std::initializer_list<label> list);

        //- Construct from FixedList of four point labels
        inline tetCell(const FixedList<label, 4>& list);

        //- Copy construct from a subset of point labels
        inline tetCell
        (
            const labelUList& list,
            const FixedList<label, 4>& indices
        );

        //- Construct from Istream
        inline tetCell(Istream& is);


    // Member Functions

        // Access

            //- Return i-th face
            inline Foam::triFace face(const label facei) const;

            //- Return first face adjacent to the given edge
            inline label edgeFace(const label edgei) const;

            //- Return face adjacent to the given face sharing the same edge
            inline label edgeAdjacentFace
            (
                const label edgei,
                const label facei
            ) const;

            //- Return i-th edge from tet.
            inline Foam::edge tetEdge(const label edgei) const;

            //- Return i-th edge. Identical to tetEdge but with generic name
            inline Foam::edge edge(const label edgei) const;

            //- Return i-th edge reversed
            inline Foam::edge reverseEdge(const label edgei) const;


        // Operations

            //- Return tet shape cell
            cellShape tetCellShape() const;

            //- Return the tetrahedron
            inline tetPointRef tet(const UList<point>& points) const;
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for tetCell
template<> struct is_contiguous<tetCell> : std::true_type {};

//- Contiguous label data for tetCell
template<> struct is_contiguous_label<tetCell> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tetCellI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
