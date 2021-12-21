/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2021 PCOpt/NTUA
    Copyright (C) 2013-2021 FOSS GP
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
    Foam::incompressible::sensitivitySurface

Description
    Calculation of adjoint based sensitivities at wall faces

    Reference:
    \verbatim
        The computation of the sensitivity derivatives  follows the (E-)SI
        formulation of

            Kavvadias, I. S., Papoutsis-Kiachagias, E. M., 
            Giannakoglou, K. C. (2015).
            On the proper treatment of grid sensitivities in continuous adjoint
            methods for shape optimization.
            Journal of computational physics, 301, 1-18.
            https://doi.org/10.1016/j.jcp.2015.08.012

        whereas their smoothing based on the computation of the 'Sobolev
        gradient' is derived from
    
            Vassberg J. C., Jameson A. (2006).
            Aerodynamic Shape Optimization Part I: Theoretical Background.
            VKI Lecture Series, Introduction to Optimization and
            Multidisciplinary Design, Brussels, Belgium, 8 March, 2006.
    \endverbatim

SourceFiles
    sensitivitySurface.C

\*---------------------------------------------------------------------------*/

#ifndef sensitivitySurfaceIncompressible_H
#define sensitivitySurfaceIncompressible_H

#include "adjointSensitivityIncompressible.H"
#include "shapeSensitivitiesBase.H"
#include "adjointEikonalSolverIncompressible.H"
#include "adjointMeshMovementSolverIncompressible.H"
#include "deltaBoundary.H"
#include "faMesh.H"

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
    public adjointSensitivity,
    public shapeSensitivitiesBase
{
protected:

    // Protected data


        //- Include surface area in sens computation
        bool includeSurfaceArea_;

        //- Include the adjoint pressure term in sens computation
        bool includePressureTerm_;

        //- Include the term containing the grad of the stress at the boundary
        bool includeGradStressTerm_;

        //- Include the transpose part of the adjoint stresses
        bool includeTransposeStresses_;

        //- Use snGrad in the transpose part of the adjoint stresses
        bool useSnGradInTranposeStresses_;

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

        //- Smooth sensitivity derivatives based on the computation of the 
        //- 'Sobolev gradient'
        bool smoothSensitivities_;

        autoPtr<adjointEikonalSolver> eikonalSolver_;

        autoPtr<adjointMeshMovementSolver> meshMovementSolver_;

        // Export face normal and face centre for use by external users
        autoPtr<volVectorField> nfOnPatchPtr_;
        autoPtr<volVectorField> SfOnPatchPtr_;
        autoPtr<volVectorField> CfOnPatchPtr_;


    // Protected Member Functions

        //- Add sensitivities from dSd/db and dnf/db computed at points and
        //- mapped to faces
        void addGeometricSens();

        //- Set suffix name for sensitivity fields
        void setSuffixName();

        //- Smooth sensitivity derivatives based on the computation of the 
        //- 'Sobolev gradient'
        void smoothSensitivities();

        //- Compute the physical smoothing radius based on the average boundary
        //- face 'length'
        scalar computeRadius(const faMesh& aMesh);


private:

    // Private Member Functions

        //- No copy construct
        sensitivitySurface(const sensitivitySurface&) = delete;

        //- No copy assignment
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
            objectiveManager& objectiveManager
        );


    //- Destructor
    virtual ~sensitivitySurface() = default;


    // Member Functions

        //- Read controls and update solver pointers if necessary
        void read();

       //- Read dict if changed
       virtual bool readDict(const dictionary& dict);

       //- Compute the number of faces on sensitivityPatchIDs_
       void computeDerivativesSize();

       //- Accumulate sensitivity integrands
       virtual void accumulateIntegrand(const scalar dt);

       //- Assemble sensitivities
       virtual void assembleSensitivities();

       //- Zero sensitivity fields and their constituents
       virtual void clearSensitivities();

       //- Get adjoint eikonal solver
       autoPtr<adjointEikonalSolver>& getAdjointEikonalSolver();

       //- Write sensitivity maps
       virtual void write(const word& baseName = word::null);

       // Inline getters and setters

           //- Get access to the includeObjective bool
           inline bool getIncludeObjective() const;

           //- Get access to the includeSurfaceArea bool
           inline bool getIncludeSurfaceArea() const;

           //- Set includeObjective bool
           inline void setIncludeObjective(const bool includeObjective);

           //- Set includeSurfaceArea bool
           inline void setIncludeSurfaceArea(const bool includeSurfaceArea);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sensitivitySurfaceIncompressibleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
