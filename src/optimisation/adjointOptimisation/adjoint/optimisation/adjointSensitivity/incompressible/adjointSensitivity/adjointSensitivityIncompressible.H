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
    Foam::incompressible::adjointSensitivity

Description
    Abstract base class for adjoint-based sensitivities in incompressible flows

    Reference:
    \verbatim
        For the FI and ESI formulations
            Kavvadias, I., Papoutsis-Kiachagias, E., & Giannakoglou, K. (2015).
            On the proper treatment of grid sensitivities in continuous adjoint
            methods for shape optimization.
            Journal of Computational Physics, 301, 1–18.
            http://doi.org/10.1016/j.jcp.2015.08.012

        For the SI formulation
            Papoutsis-Kiachagias, E. M., & Giannakoglou, K. C. (2014).
            Continuous Adjoint Methods for Turbulent Flows, Applied to Shape
            and Topology Optimization: Industrial Applications.
            Archives of Computational Methods in Engineering, 23(2), 255–299.
            http://doi.org/10.1007/s11831-014-9141-9
    \endverbatim

SourceFiles
    adjointSensitivity.C

\*---------------------------------------------------------------------------*/

#ifndef adjointSensitivityIncompressible_H
#define adjointSensitivityIncompressible_H

#include "sensitivity.H"
#include "incompressibleVars.H"
#include "incompressibleAdjointVars.H"
#include "fvOptionAdjointList.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace incompressible
{

/*---------------------------------------------------------------------------*\
                     Class adjointSensitivity Declaration
\*---------------------------------------------------------------------------*/

class adjointSensitivity
:
    public sensitivity
{
protected:

    // Protected data

        incompressibleVars& primalVars_;
        incompressibleAdjointVars& adjointVars_;
        objectiveManager& objectiveManager_;
        fv::optionAdjointList& fvOptionsAdjoint_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        adjointSensitivity(const adjointSensitivity&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const adjointSensitivity&) = delete;


public:

    //- Runtime type information
    TypeName("adjointSensitivity");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            adjointSensitivity,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict,
                incompressibleVars& primalVars,
                incompressibleAdjointVars& adjointVars,
                objectiveManager& objectiveManager,
                fv::optionAdjointList& fvOptionsAdjoint
            ),
            (
                mesh,
                dict,
                primalVars,
                adjointVars,
                objectiveManager,
                fvOptionsAdjoint
            )
        );


    // Constructors

        //- Construct from components
        adjointSensitivity
        (
            const fvMesh& mesh,
            const dictionary& dict,
            incompressibleVars& primalVars,
            incompressibleAdjointVars& adjointVars,
            objectiveManager& objectiveManager,
            fv::optionAdjointList& fvOptionsAdjoint
        );

    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<adjointSensitivity> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            incompressibleVars& primalVars,
            incompressibleAdjointVars& adjointVars,
            objectiveManager& objectiveManager,
            fv::optionAdjointList& fvOptionsAdjoint
        );


    //- Destructor
    virtual ~adjointSensitivity() = default;


    // Member Functions

       //- Calculates and returns sensitivity fields.
       //  Used with optimisation libraries
       virtual const scalarField& calculateSensitivities() = 0;

       //- Write sensitivity fields.
       //  If valid, copies boundaryFields to volFields and writes them.
       //  Virtual to be reimplemented by control points-based methods
       //  (Bezier, RBF) which do not need to write fields
       virtual void write(const word& baseName = word::null);

       //- Compute the volTensorField multiplying grad(dxdb) for
       //- the volume-based approach to compute shape sensitivity derivatives
       tmp<volTensorField> computeGradDxDbMultiplier();

       //- Compute source term for adjoint mesh movement equation
       tmp<volVectorField> adjointMeshMovementSource();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
