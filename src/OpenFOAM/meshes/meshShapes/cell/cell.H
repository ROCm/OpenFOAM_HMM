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

Class
    Foam::cell

Description
    A cell is defined as a list of faces with extra functionality

SourceFiles
    cellI.H
    cell.C
    oppositeCellFace.C

\*---------------------------------------------------------------------------*/

#ifndef cell_H
#define cell_H

#include "faceList.H"
#include "oppositeFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class cell;
bool operator==(const cell& a, const cell& b);
inline bool operator!=(const cell& a, const cell& b);


/*---------------------------------------------------------------------------*\
                           Class cell Declaration
\*---------------------------------------------------------------------------*/

class cell
:
    public labelList
{

public:

    // Static data members

        static const char* const typeName;


    // Constructors

        //- Construct null
        inline cell();

        //- Construct given size
        explicit inline cell(const label sz);

        //- Construct from list of labels
        explicit inline cell(const labelUList& lst);

        //- Construct by transferring list of labels
        explicit inline cell(const Xfer<labelList>& lst);

        //- Construct from Istream
        inline cell(Istream& is);


    // Member Functions

        //- Return number of faces
        inline label nFaces() const;

        //- Return unordered list of cell vertices given the list of faces
        labelList labels(const faceUList& f) const;

        //- Return the cell vertices given the list of faces and mesh points
        pointField points
        (
            const faceUList& f,
            const UList<point>& meshPoints
        ) const;

        //- Return cell edges
        edgeList edges(const faceUList& f) const;

        //- Return index of opposite face
        label opposingFaceLabel
        (
            const label masterFaceLabel,
            const faceUList& meshFaces
        ) const;

        //- Return opposite face oriented the same way as the master face
        oppositeFace opposingFace
        (
            const label masterFaceLabel,
            const faceUList& meshFaces
        ) const;


        // Instructions for centre and mag:
        // When one wants to access the cell centre and magnitude, the
        // functionality on the mesh level should be used in preference to the
        // functions provided here. They do not rely to the functionality
        // implemented here, provide additional checking and are more efficient.
        // The cell::centre and cell::mag functions may be removed in the
        // future.

        //- Returns cell centre
        point centre(const UList<point>& p, const faceUList& f) const;

        //- Returns cell volume
        scalar mag(const UList<point>& p, const faceUList& f) const;


    // Friend Operators

        friend bool operator==(const cell& a, const cell& b);
        friend bool operator!=(const cell& a, const cell& b);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
