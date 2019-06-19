/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2007-2019 PCOpt/NTUA
                            | Copyright (C) 2013-2019 FOSS GP
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
    Foam::incompressible::sensitivitySurface

Description
    Calculation of adjoint based sensitivities at wall faces

SourceFiles
    sensitivitySurface.C

\*---------------------------------------------------------------------------*/

#ifndef sensitivitySurfaceIncompressible_H
#define sensitivitySurfaceIncompressible_H

#include "adjointSensitivityIncompressible.H"
#include "adjointEikonalSolverIncompressible.H"
#include "adjointMeshMovementSolverIncompressible.H"
#include "deltaBoundary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace incompressible
{

/*---------------------------------------------------------------------------*\
                     Class sensitivitySurface Declaration
\*---------------------------------------------------------------------------*/

class sensitivitySurface
:
    public adjointSensitivity
{
protected:

    // Protected data

        //- Scalar normal sens
        scalarField derivatives_;

        //- Include surface area in sens computation
        bool includeSurfaceArea_;

        //- Include the adjoint pressure term in sens computation
        bool includePressureTerm_;

        //- Include the term containing the grad of the stress at the boundary
        bool includeGradStressTerm_;

        //- Include the transpose part of the adjoint stresses
        bool includeTransposeStresses_;

        //- Include the term from the deviatoric part of the stresses
        bool includeDivTerm_;

        //- Include distance variation in sens computation
        bool includeDistance_;

        //- Include mesh movement variation in sens computation
        bool includeMeshMovement_;

        //- Include terms directly emerging from the objective function
        bool includeObjective_;

        //- Write geometric info for use by external programs
        bool writeGeometricInfo_;

        autoPtr<adjointEikonalSolver> eikonalSolver_;

        autoPtr<adjointMeshMovementSolver> meshMovementSolver_;

        // Export face normal and face centre for use by external users
        autoPtr<volVectorField> nfOnPatchPtr_;
        autoPtr<volVectorField> SfOnPatchPtr_;
        autoPtr<volVectorField> CfOnPatchPtr_;


    // Protected Member Functions

        //- Read controls and update solver pointers if necessary
        void read();

        //- Add sensitivities from dSd/db and dnf/db computed at points and
        //- mapped to faces
        void addGeometricSens();


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        sensitivitySurface(const sensitivitySurface&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const sensitivitySurface&) = delete;


public:

    //- Runtime type information
    TypeName("surface");


    // Constructors

        //- Construct from components
        sensitivitySurface
        (
            const fvMesh& mesh,
            const dictionary& dict,
            incompressibleVars& primalVars,
            incompressibleAdjointVars& adjointVars,
            objectiveManager& objectiveManager,
            fv::optionAdjointList& fvOptionsAdjoint
        );


    //- Destructor
    virtual ~sensitivitySurface() = default;


    // Member Functions

       //- Read dict if changed
       virtual bool readDict(const dictionary& dict);

       //- Compute the number of faces on sensitivityPatchIDs_
       void computeDerivativesSize();

       //- Calculates sensitivities at wall surface points
       const scalarField& calculateSensitivities();

       //- Get adjoint eikonal solver
       autoPtr<adjointEikonalSolver>& getAdjointEikonalSolver();

       virtual void write(const word& baseName = word::null);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
