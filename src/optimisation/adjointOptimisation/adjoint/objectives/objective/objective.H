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
    Foam::objective

Description
    Abstract base class for objective functions. No point in making this
    runTime selectable since its children will have different constructors.

SourceFiles
    objective.C

\*---------------------------------------------------------------------------*/

#ifndef objective_H
#define objective_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "OFstream.H"
#include "boundaryFieldsFwd.H"
#include "solverControl.H"
#include "objectiveFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class objective Declaration
\*---------------------------------------------------------------------------*/

class objective
{
protected:

    // Protected data

        const fvMesh& mesh_;
        dictionary dict_;
        const word adjointSolverName_;
        const word primalSolverName_;
        const word objectiveName_;
        bool computeMeanFields_;

        // Objective function value and weight
        scalar J_;
        scalar JMean_;  //average value
        scalar weight_;

        // Contribution to field sensitivity derivatives
        // Topology optimisation
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        autoPtr<volScalarField> dJdbPtr_;

        // Contribution to surface sensitivity derivatives
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        //- Term from material derivative
        autoPtr<boundaryVectorField> bdJdbPtr_;

        //- Term multiplying delta(n dS)/delta b
        autoPtr<boundaryVectorField> bdSdbMultPtr_;

        //- Term multiplying delta(n)/delta b
        autoPtr<boundaryVectorField> bdndbMultPtr_;

        //- Term multiplying delta(x)/delta b at the boundary
        autoPtr<boundaryVectorField> bdxdbMultPtr_;

        //- Term multiplying delta(x)/delta b at the boundary
        //- for objectives that directly depend on x, e.g. moment
        //- Needed in both FI and SI computations
        autoPtr<boundaryVectorField> bdxdbDirectMultPtr_;

        //- Contribution located in specific parts of a patch.
        //- Mainly intended for patch boundary edges contributions, e.g.
        //- NURBS surfaces G1 continuity
        autoPtr<vectorField3> bEdgeContribution_;

        //- For use with discrete-like sensitivities
        autoPtr<boundaryTensorField> bdJdStressPtr_;

        // Contribution to volume-based sensitivities from volume-based
        // objective functions
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        //- Multiplier of d(Volume)/db
        autoPtr<volScalarField> divDxDbMultPtr_;

        //- Emerging from volume objectives that include spatial derivatives
        autoPtr<volTensorField> gradDxDbMultPtr_;

        //- Output file variables
        fileName objFunctionFolder_;

        //- File to keep the objective values after the end of the primal solver
        mutable autoPtr<OFstream> objFunctionFilePtr_;

        //- File to keep the objective values at each iteration of the primal
        //- solver
        mutable autoPtr<OFstream> instantValueFilePtr_;

        //- File to keep the average objective values after the end of the
        //- primal solver
        mutable autoPtr<OFstream> meanValueFilePtr_;


    // Protected Member Functions

        //- Return objective dictionary
        const dictionary& dict() const;

        //- Set the output file ptr
        void setObjectiveFilePtr() const;

        //- Set the output file ptr for the instantaneous value
        void setInstantValueFilePtr() const;

        //- Set the output file ptr for the mean value
        void setMeanValueFilePtr() const;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        objective(const objective&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const objective&) = delete;

        //- Make objective Function Folder
        void makeFolder();


public:

    //- Runtime type information
    TypeName("objective");


    // Declare run-time constructor selection table

        declareRunTimeNewSelectionTable
        (
            autoPtr,
            objective,
            objective,
            (
                const fvMesh& mesh,
                const dictionary& dict,
                const word& adjointSolverName,
                const word& primalSolverName
            ),
            (mesh, dict, adjointSolverName, primalSolverName)
        );


    // Constructors

        //- Construct from components
        objective
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& adjointSolverName,
            const word& primalSolverName
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<objective> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& objectiveType,
            const word& adjointSolverName,
            const word& primalSolverName
        );


    //- Destructor
    virtual ~objective() = default;


    // Member Functions

        virtual bool readDict(const dictionary& dict);

        //- Return the objective function value
        virtual scalar J() = 0;

        //- Accumulate contribution for the mean objective value
        void accumulateJMean(solverControl& solverControl);

        //- Return the objective function weight
        scalar weight() const;

        //- Contribution to field sensitivities
        const volScalarField& dJdb();

        //- Contribution to surface sensitivities for a specific patch
        const fvPatchVectorField& boundarydJdb(const label);

        //- Multiplier of delta(n dS)/delta b
        const fvPatchVectorField& dSdbMultiplier(const label);

        //- Multiplier of delta(n dS)/delta b
        const fvPatchVectorField& dndbMultiplier(const label);

        //- Multiplier of delta(x)/delta b
        const fvPatchVectorField& dxdbMultiplier(const label);

        //- Multiplier of delta(x)/delta b
        const fvPatchVectorField& dxdbDirectMultiplier(const label);

        //- Multiplier located at patch boundary edges
        const vectorField& boundaryEdgeMultiplier
        (
            const label patchI,
            const label edgeI
        );

        //- Objective partial deriv wrt stress tensor
        const fvPatchTensorField& boundarydJdStress(const label);

        //- Contribution to surface sensitivities for all patches
        const boundaryVectorField& boundarydJdb();

        //- Multiplier of delta(n dS)/delta b for all patches
        const boundaryVectorField& dSdbMultiplier();

        //- Multiplier of delta(n dS)/delta b for all patches
        const boundaryVectorField& dndbMultiplier();

        //- Multiplier of delta(x)/delta b for all patches
        const boundaryVectorField& dxdbMultiplier();

        //- Multiplier of delta(x)/delta b for all patches
        const boundaryVectorField& dxdbDirectMultiplier();

        //- Multiplier located at patch boundary edges
        const vectorField3& boundaryEdgeMultiplier();

        //- Objective partial deriv wrt stress tensor
        const boundaryTensorField& boundarydJdStress();

        //- Multiplier of grad( delta(x)/delta b) for volume-based sensitivities
        const volScalarField& divDxDbMultiplier();

        //- Multiplier of grad( delta(x)/delta b) for volume-based sensitivities
        const volTensorField& gradDxDbMultiplier();

        //- Update objective function derivatives
        virtual void update() = 0;

        //- Update normalization factors, for objectives in
        //- which the factor is not known a priori
        virtual void updateNormalizationFactor();

        //- Update objective function derivative term
        virtual void update_boundarydJdb()
        {}

        //- Update d (normal dS) / db multiplier. Surface-based sensitivity term
        virtual void update_dSdbMultiplier()
        {}

        //- Update d (normal) / db multiplier. Surface-based sensitivity term
        virtual void update_dndbMultiplier()
        {}

        //- Update d (x) / db multiplier. Surface-based sensitivity term
        virtual void update_dxdbMultiplier()
        {}

        //- Update d (x) / db multiplier. Surface and volume-based sensitivity
        //- term
        virtual void update_dxdbDirectMultiplier()
        {}

        //- Update boundary edge contributions
        virtual void update_boundaryEdgeContribution()
        {}

        //- Update dJ/dStress field
        virtual void update_dJdStressMultiplier()
        {}

        //- Update div( dx/db multiplier). Volume-based sensitivity term
        virtual void update_divDxDbMultiplier()
        {}

        //- Update grad( dx/db multiplier). Volume-based sensitivity term
        virtual void update_gradDxDbMultiplier()
        {}

        //- Write objective function history
        virtual void write() const;

        //- Write objective function history at each primal solver iteration
        virtual void writeInstantaneousValue() const;

        //- Write mean objective function history
        virtual void writeMeanValue() const;

        // Inline functions for checking whether pointers are set or not
        inline const word& objectiveName() const;
        inline bool hasdJdb();
        inline bool hasBoundarydJdb();
        inline bool hasdSdbMult();
        inline bool hasdndbMult();
        inline bool hasdxdbMult();
        inline bool hasdxdbDirectMult();
        inline bool hasBoundaryEdgeContribution();
        inline bool hasBoundarydJdStress();
        inline bool hasDivDxDbMult();
        inline bool hasGradDxDbMult();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "objectiveI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
