/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::oversetLduInterface

Description
    An abstract base class for overset coupled interfaces

SourceFiles
    oversetLduInterface.C

\*---------------------------------------------------------------------------*/

#ifndef oversetLduInterface_H
#define oversetLduInterface_H

#include "mapDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class oversetLduInterface Declaration
\*---------------------------------------------------------------------------*/

class oversetLduInterface
{

public:

    //- Runtime type information
    TypeName("oversetLduInterface");


    // Constructors

        //- Construct null
        oversetLduInterface()
        {}


    //- Destructor
    virtual ~oversetLduInterface();


    // Member Functions

        // Access

            //- Name of interface (for debugging)
            virtual const word& name() const = 0;

            //- Am I the master interface
            virtual bool master() const = 0;

            //- GAMG restriction (fine-to-coarse)
            virtual const labelList& restrictMap() const = 0;

            //- Donor stencil
            virtual const labelListList& stencil() const = 0;

            //- Map for obtaining data in stencil order
            virtual const mapDistribute& cellInterpolationMap() const = 0;

            //- Weights in stencil order
            virtual const List<scalarList>&
            cellInterpolationWeights() const = 0;

            //- Normalisation of matrix; for explicit contributions
            virtual const scalarField& normalisation() const = 0;

            //- Acceptor cells
            virtual const labelList& interpolationCells() const = 0;

            //- Underrelaxation for acceptor cells
            virtual const scalarList& cellInterpolationWeight() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
