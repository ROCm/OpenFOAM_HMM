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
    Foam::incompressible::RASModelVariables

Description
    Abstract base class for objective functions. No point in making this
    runTime selectable since its childs will have different constructors.

SourceFiles
    RASModelVariables.C

\*---------------------------------------------------------------------------*/

#ifndef RASModelVariables_H
#define RASModelVariables_H

#include "solverControl.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

/*---------------------------------------------------------------------------*\
                      Class RASModelVariables Declaration
\*---------------------------------------------------------------------------*/

class RASModelVariables
{

public:

    // Public typedefs

        typedef autoPtr<tmp<volScalarField>> autoTmp;

protected:

    // Protected data


        const fvMesh& mesh_;
        const solverControl& solverControl_;

        // autoPtrs delete the memory on destruction
        // Can cause memory mishandling issues in this case
        // Use regular ptrs instead
        bool hasTMVar1_;
        bool hasTMVar2_;
        bool hasNut_;
        bool hasDist_;
        autoTmp TMVar1Ptr_;
        autoTmp TMVar2Ptr_;
        autoTmp nutPtr_;
        autoTmp dPtr_;
        // Base names of the turbulent fields
        word TMVar1BaseName_;
        word TMVar2BaseName_;
        word nutBaseName_;

        // conditionally store initial values
        // For finite differences and optimisation runs
        autoPtr<volScalarField> TMVar1InitPtr_;
        autoPtr<volScalarField> TMVar2InitPtr_;
        autoPtr<volScalarField> nutInitPtr_;

        // conditionally store mean values
        autoPtr<volScalarField> TMVar1MeanPtr_;
        autoPtr<volScalarField> TMVar2MeanPtr_;
        autoPtr<volScalarField> nutMeanPtr_;


    // Protected Member functions

        void allocateInitValues();
        void allocateMeanFields();
        autoTmp cloneAutoTmp(const autoTmp& source) const;
        void copyAndRename(volScalarField& f1, volScalarField& f2);


private:

    // Private Member Functions


        //- No copy assignment
        void operator=(const RASModelVariables&) = delete;


public:


    //- Runtime type information
    TypeName("RASModelVariables");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            RASModelVariables,
            dictionary,
            (
                const fvMesh& mesh,
                const solverControl& SolverControl
            ),
            (mesh, SolverControl)
        );


    // Constructors

        //- Construct from components
        RASModelVariables
        (
            const fvMesh& mesh,
            const solverControl& SolverControl
        );

        //- Copy constructor
        //  Will allocate new fields (instead of referencing the ones in the
        //  turbulence model), so cannot be used directly to access the fields
        //  of the turbulence model.
        //  Mainly used for checkpointing in unsteady adjoint
        RASModelVariables
        (
            const RASModelVariables& rmv
        );

        //- Clone
        //  Will allocate new fields (instead of referencing the ones in the
        //  turbulence model), so cannot be used directly to access the fields
        //  of the turbulence model.
        //  Mainly used for checkpointing in unsteady adjoint
        autoPtr<RASModelVariables> clone() const;


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<RASModelVariables> New
        (
            const fvMesh& mesh,
            const solverControl& SolverControl
        );


    // Destructor

        // Destructor does nothing on base since depending on the case new
        // fields might be allocated
        // MUST be overloaded in inherited classes
        virtual ~RASModelVariables() = default;


    // Member Functions


        //- Bools to idenify which turbulent fields are present
        bool hasTMVar1() const;
        bool hasTMVar2() const;
        bool hasNut() const;
        bool hasDist() const;

        //- Turbulence field names
        const word& TMVar1BaseName() const;
        const word& TMVar2BaseName() const;
        const word& nutBaseName() const;

        //- Return references to turbulence fields
        //  will return the mean field if it exists, otherwise the
        //  instantaneous one
        const volScalarField& TMVar1() const;
              volScalarField& TMVar1();
        const volScalarField& TMVar2() const;
              volScalarField& TMVar2();
        const volScalarField& nutRef() const;
              volScalarField& nutRef();
        const volScalarField& d() const;
              volScalarField& d();

        //- return references to instantaneous turbulence fields
        const volScalarField& TMVar1Inst() const;
              volScalarField& TMVar1Inst();
        const volScalarField& TMVar2Inst() const;
              volScalarField& TMVar2Inst();
        const volScalarField& nutRefInst() const;
              volScalarField& nutRefInst();

        //- Return nut Jacobian wrt the TM vars
        virtual tmp<volScalarField> nutJacobianVar1
        (
            const singlePhaseTransportModel& laminarTransport
        ) const;
        virtual tmp<volScalarField> nutJacobianVar2
        (
            const singlePhaseTransportModel& laminarTransport
        ) const;

        //- Restore turbulent fields to their initial values
        void restoreInitValues();

        //- Reset mean fields to zero
        void resetMeanFields();

        //- Compute mean fields on the fly
        void computeMeanFields();

        //- Return stress tensor based on the mean flow variables
        tmp<volSymmTensorField> devReff
        (
            const singlePhaseTransportModel& laminarTransport,
            const volVectorField& U
        ) const;

        //- correct bounday conditions of turbulent fields
        virtual void correctBoundaryConditions
        (
            const incompressible::turbulenceModel& turbulence
        );

        //- Transfer turbulence fields from an another object
        //  Copies values since the ownership of the original fields is held
        //  by the turbulence model
        virtual void transfer(RASModelVariables& rmv);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
