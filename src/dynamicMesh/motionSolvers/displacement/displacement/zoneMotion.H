/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::zoneMotion

Description


SourceFiles
    zoneMotion.C

\*---------------------------------------------------------------------------*/

#ifndef zoneMotion_H
#define zoneMotion_H

#include "labelList.H"
#include "dictionary.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class zoneMotion Declaration
\*---------------------------------------------------------------------------*/

class zoneMotion
{
     // Private Data

        //- Points to move when cell zone is supplied
        labelList pointIDs_;

        //- Flag to indicate whether all cells should move
        bool moveAllCells_;


    // Private Member Functions

        //- No copy construct
        zoneMotion(const zoneMotion&) = delete;

        //- No copy assignment
        void operator=(const zoneMotion&) = delete;

public:

    // Member functions

        //- Return pointsID
        const labelList& pointIDs() const;

        //- Return flag
        bool moveAllCells() const;

    // Constructors

        //- Construct from dictionary
        zoneMotion
        (
            const dictionary&,
            const polyMesh&
        );


    //- Destructor
    virtual ~zoneMotion() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
