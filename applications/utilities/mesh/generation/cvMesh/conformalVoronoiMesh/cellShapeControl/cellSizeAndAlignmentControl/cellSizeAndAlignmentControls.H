/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::cellSizeAndAlignmentControls

Description

SourceFiles
    cellSizeAndAlignmentControls.C

\*---------------------------------------------------------------------------*/

#ifndef cellSizeAndAlignmentControls_H
#define cellSizeAndAlignmentControls_H

#include "dictionary.H"
#include "cellShapeControlMesh.H"
#include "cellSizeAndAlignmentControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                     Class cellSizeAndAlignmentControls Declaration
\*---------------------------------------------------------------------------*/

class cellSizeAndAlignmentControls
{
    // Private data

        const dictionary& shapeControlDict_;

        const conformationSurfaces& allGeometry_;

        PtrList<cellSizeAndAlignmentControl> controlFunctions_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        cellSizeAndAlignmentControls(const cellSizeAndAlignmentControls&);

        //- Disallow default bitwise assignment
        void operator=(const cellSizeAndAlignmentControls&);


public:

    // Constructors

        //- Construct from dictionary
        cellSizeAndAlignmentControls
        (
            const Time& runTime,
            const dictionary& shapeControlDict,
            const conformationSurfaces& allGeometry
        );


    //- Destructor
    virtual ~cellSizeAndAlignmentControls();


    // Member Functions

        // Access

            const PtrList<cellSizeAndAlignmentControl>& controlFunctions() const
            {
                return controlFunctions_;
            }


        // Query
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
