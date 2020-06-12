/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::objectiveIncompressible

Description
    Abstract base class for objective functions in incompressible flows

SourceFiles
    objectiveIncompressible.C

\*---------------------------------------------------------------------------*/

#ifndef objectiveIncompressible_H
#define objectiveIncompressible_H

#include "objective.H"
#include "incompressibleVars.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class objectiveIncompressible Declaration
\*---------------------------------------------------------------------------*/

class objectiveIncompressible
:
    public objective
{
protected:

    // Protected data

        const incompressibleVars& vars_;

        // Contribution to field adjoint equations
        // v,p,T and turbulence model variables
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        autoPtr<volVectorField> dJdvPtr_;
        autoPtr<volScalarField> dJdpPtr_;
        autoPtr<volScalarField> dJdTPtr_;

        //- First  turbulence model variable
        autoPtr<volScalarField> dJdTMvar1Ptr_;

        //- Second turbulence model variable
        autoPtr<volScalarField> dJdTMvar2Ptr_;

        // Contribution to adjoint boundary conditions
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        autoPtr<boundaryVectorField> bdJdvPtr_;

        //- Adjoint outlet pressure
        autoPtr<boundaryScalarField> bdJdvnPtr_;

        //- Adjoint outlet velocity
        autoPtr<boundaryVectorField> bdJdvtPtr_;

        //- Adjoint (intlet,wall) velocity
        autoPtr<boundaryVectorField> bdJdpPtr_;

        //- Adjoint outlet temperature
        autoPtr<boundaryScalarField> bdJdTPtr_;

        //- Adjoint outlet turbulence model var 1
        autoPtr<boundaryScalarField> bdJdTMvar1Ptr_;

        //- Adjoint outlet turbulence model var 2
        autoPtr<boundaryScalarField> bdJdTMvar2Ptr_;


private:

    // Private Member Functions

        //- No copy construct
        objectiveIncompressible(const objectiveIncompressible&) = delete;

        //- No copy assignment
        void operator=(const objectiveIncompressible&) = delete;


public:

    //- Runtime type information
    TypeName("incompressible");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            objectiveIncompressible,
            dictionary,
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
        objectiveIncompressible
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& adjointSolverName,
            const word& primalSolverName
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<objectiveIncompressible> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& adjointSolverName,
            const word& primalSolverName
        );


    //- Destructor
    virtual ~objectiveIncompressible() = default;


    // Member Functions

        //- Return the objective function value
        virtual scalar J() = 0;

        //- Normalize all fields allocated by the objective
        virtual void doNormalization();

        //- Contribution to field adjoint momentum eqs
        const volVectorField& dJdv();

        //- Contribution to field adjoint continuity eq
        const volScalarField& dJdp();

        //- Contribution to field adjoint energy eq
        const volScalarField& dJdT();

        //- Contribution to field adjoint turbulence model variable 1
        const volScalarField& dJdTMvar1();

        //- Contribution to field adjoint turbulence model variable 2
        const volScalarField& dJdTMvar2();

        //- Objective partial deriv wrt velocity for a specific patch
        const fvPatchVectorField& boundarydJdv(const label);

        //- Objective partial deriv wrt normal velocity for a specific patch
        const fvPatchScalarField& boundarydJdvn(const label);

        //- Objective partial deriv wrt tangent velocity for a specific patch
        const fvPatchVectorField& boundarydJdvt(const label);

        //- Objective partial deriv wrt pressure (times normal) for a specific
        //- patch
        const fvPatchVectorField& boundarydJdp(const label);

        //- Objective partial deriv wrt temperature for a specific patch
        const fvPatchScalarField& boundarydJdT(const label);

        //- Objective partial deriv wrt turbulence model var 1 for a specific
        //- patch
        const fvPatchScalarField& boundarydJdTMvar1(const label);

        //- Objective partial deriv wrt turbulence model var 2 for a specific
        //- patch
        const fvPatchScalarField& boundarydJdTMvar2(const label);

        //- Objective partial deriv wrt velocity for all patches
        const boundaryVectorField& boundarydJdv();

        //- Objective partial deriv wrt normal velocity for all patches
        const boundaryScalarField& boundarydJdvn();

        //- Objective partial deriv wrt tangent velocity for all patches
        const boundaryVectorField& boundarydJdvt();

        //- Objective partial deriv wrt pressure (times normal) for all patches
        const boundaryVectorField& boundarydJdp();

        //- Objective partial deriv wrt temperature for all patches
        const boundaryScalarField& boundarydJdT();

        //- Objective partial deriv wrt turbulence model var 1 for all patches
        const boundaryScalarField& boundarydJdTMvar1();

        //- Objective partial deriv wrt turbulence model var 2 for all patches
        const boundaryScalarField& boundarydJdTMvar2();

        //- Update objective function derivatives
        virtual void update();

        //- Update objective function derivatives
        virtual void nullify();

        //- Update vol and boundary fields and derivatives
        //  Do nothing in the base. The relevant ones should be overwritten
        //  in the child objective functions
        virtual void update_dJdv()
        {}

        virtual void update_dJdp()
        {}

        virtual void update_dJdT()
        {}

        virtual void update_dJdTMvar1()
        {}

        virtual void update_dJdTMvar2()
        {}

        virtual void update_dJdb()
        {}

        virtual void update_divDxDbMultiplier()
        {}

        virtual void update_gradDxDbMultiplier()
        {}

        virtual void update_boundarydJdv()
        {}

        virtual void update_boundarydJdvn()
        {}

        virtual void update_boundarydJdvt()
        {}

        virtual void update_boundarydJdp()
        {}

        virtual void update_boundarydJdT()
        {}

        virtual void update_boundarydJdTMvar1()
        {}

        virtual void update_boundarydJdTMvar2()
        {}

        virtual void update_boundarydJdb()
        {}

        virtual void update_dSdbMultiplier()
        {}

        virtual void update_dndbMultiplier()
        {}

        virtual void update_dxdbMultiplier()
        {}

        virtual void update_dxdbDirectMultiplier()
        {}

        //- Some objectives need to store some auxiliary values.
        //- If averaging is enabled, update these mean values here.
        //  By convention, the mean values (eg mean drag) refer to these flow
        //  values computed using the mean fields, rather than averaging the
        //  values themselves
        virtual void update_meanValues()
        {}

        //- Write objective function history
        virtual bool write(const bool valid = true) const;

        //- Inline functions for checking whether pointers are set or not
        inline bool hasdJdv() const;
        inline bool hasdJdp() const;
        inline bool hasdJdT() const;
        inline bool hasdJdTMVar1() const;
        inline bool hasdJdTMVar2() const;
        inline bool hasBoundarydJdv() const;
        inline bool hasBoundarydJdvn() const;
        inline bool hasBoundarydJdvt() const;
        inline bool hasBoundarydJdp() const;
        inline bool hasBoundarydJdT() const;
        inline bool hasBoundarydJdTMVar1() const;
        inline bool hasBoundarydJdTMVar2() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "objectiveIncompressibleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
