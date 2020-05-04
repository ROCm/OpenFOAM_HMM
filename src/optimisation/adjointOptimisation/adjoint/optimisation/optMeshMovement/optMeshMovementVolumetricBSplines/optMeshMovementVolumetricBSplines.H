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
    Foam::optMeshMovementVolumetricBSplines

Description
    Converts NURBS volume control points update to actual mesh movement.
    Internal points are also moved based on the movement of the control points

SourceFiles
    optMeshMovementVolumetricBSplines.C

\*---------------------------------------------------------------------------*/

#ifndef optMeshMovementVolumetricBSplines_H
#define optMeshMovementVolumetricBSplines_H

#include "optMeshMovement.H"
#include "volBSplinesBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class optMeshMovementVolumetricBSplines Declaration
\*---------------------------------------------------------------------------*/

class optMeshMovementVolumetricBSplines
:
    public optMeshMovement
{
protected:

    // Protected data

        //- Reference to underlaying volumetric B-Splines morpher
        volBSplinesBase& volBSplinesBase_;

        //- Backup of initial control points. Useful for line-search
        List<vectorField> cpsInit_;


    // Protected Member Functions

        vectorField controlPointMovement(const scalarField& correction);


private:

    // Private Member Functions

        //- No copy construct
        optMeshMovementVolumetricBSplines
        (
            const optMeshMovementVolumetricBSplines&
        ) = delete;

        //- No copy assignment
        void operator=(const optMeshMovementVolumetricBSplines&) = delete;


public:

    //- Runtime type information
    TypeName("volumetricBSplines");


    // Constructors

        //- Construct from components
        optMeshMovementVolumetricBSplines
        (
            fvMesh& mesh,
            const dictionary& dict,
            const labelList& patchIDs
        );


    //- Destructor
    virtual ~optMeshMovementVolumetricBSplines() = default;


    // Member Functions

       //- Calculates surface mesh movement
       void moveMesh();

       //- Store design variables and mesh, to act as the starting point of
       //- line search
       virtual void storeDesignVariables();

       //- Reset to starting point of line search
       virtual void resetDesignVariables();

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
