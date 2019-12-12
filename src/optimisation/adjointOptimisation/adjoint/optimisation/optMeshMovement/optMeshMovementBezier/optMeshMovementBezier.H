/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
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
    Foam::optMeshMovementBezier

Description
    Converts NURBS control points update to actual mesh movement

SourceFiles
    optMeshMovementBezier.C

\*---------------------------------------------------------------------------*/

#ifndef optMeshMovementBezier_H
#define optMeshMovementBezier_H

#include "optMeshMovement.H"
#include "Bezier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class optMeshMovementBezier Declaration
\*---------------------------------------------------------------------------*/

class optMeshMovementBezier
:
    public optMeshMovement
{
protected:

    // Protected data

        //- Parameterization based on NURBS curves
        Bezier Bezier_;

        //- Boundary movement due to change of NURBS control points
        pointVectorField dx_;

        //- Cumulative change of control points
        vectorField cumulativeChange_;


    // Protected Member Functions

        void computeBoundaryMovement(const scalarField& correction);


private:

    // Private Member Functions

        //- No copy construct
        optMeshMovementBezier(const optMeshMovementBezier&) = delete;

        //- No copy assignment
        void operator=(const optMeshMovementBezier&) = delete;


public:

    //- Runtime type information
    TypeName("Bezier");


    // Constructors

        //- Construct from components
        optMeshMovementBezier
        (
            fvMesh& mesh,
            const dictionary& dict,
            const labelList& patchIDs
        );


    //- Destructor
    virtual ~optMeshMovementBezier() = default;


    // Member Functions

       //- Calculates surface mesh movement
       void moveMesh();

       //- Compute eta value based on max displacement
       virtual scalar computeEta(const scalarField& correction);

       //- Return active design variables
       virtual labelList getActiveDesignVariables() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
