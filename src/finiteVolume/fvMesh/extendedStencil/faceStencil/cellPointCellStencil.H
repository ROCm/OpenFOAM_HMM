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
    Foam::cellPointCellStencil

Description

SourceFiles
    cellPointCellStencil.C

\*---------------------------------------------------------------------------*/

#ifndef cellPointCellStencil_H
#define cellPointCellStencil_H

#include "faceStencil.H"
#include "boolList.H"
#include "HashSet.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cellPointCellStencil Declaration
\*---------------------------------------------------------------------------*/

class cellPointCellStencil
:
    public faceStencil
{
    // Private Member Functions

        //- Calculates per point the neighbour data (= pointCells)
        void calcPointBoundaryData
        (
            const boolList& isValidBFace,
            const labelList& boundaryPoints,
            Map<labelList>& neiGlobal
        ) const;

        void calcCellStencil(labelListList& globalCellCells) const;


        //- Disallow default bitwise copy construct
        cellPointCellStencil(const cellPointCellStencil&);

        //- Disallow default bitwise assignment
        void operator=(const cellPointCellStencil&);


public:

    // Constructors

        //- Construct from all cells and boundary faces
        explicit cellPointCellStencil(const polyMesh&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
