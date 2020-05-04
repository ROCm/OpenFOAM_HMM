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
    Foam::updateMethod

Description
    Abstract base class for optimisation methods

SourceFiles
    updateMethod.C

\*---------------------------------------------------------------------------*/

#ifndef updateMethod_H
#define updateMethod_H

#include "runTimeSelectionTables.H"
#include "IOdictionary.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class updateMethod Declaration
\*---------------------------------------------------------------------------*/

class updateMethod
{
protected:

    // Protected data

        const fvMesh& mesh_;

        const dictionary dict_;

        //- Used to output values useful for continuation runs
        IOdictionary optMethodIODict_;

        //- Derivatives of the objective functions
        scalarField objectiveDerivatives_;

        //- Derivatives of the constraints
        PtrList<scalarField> constraintDerivatives_;

        //- Objective value
        scalar objectiveValue_;

        //- Constraint values
        scalarField cValues_;

        //- Design variables correction
        scalarField correction_;

        //- Cumulative design variables correction throughout the optimisation
        //- loop
        scalarField cumulativeCorrection_;

        //- Step multiplying the correction
        scalar eta_;

        //- Is initially set?
        bool initialEtaSet_;

        //- Folder storing the corrections to file
        //  For some optimisation methods with a very high number of
        //  design variables (e.g. topology), it doesn't make much sense
        //  to write all updates in the updateMethodDict. Hence, a
        //  separate file is used to write the corrections, in case they are
        //  needed for post-processing
        word correctionFolder_;

        //- Whether to use gSum or sum in the inner products
        bool globalSum_;

        // Scalar -- matrix multiplications
        const scalarField leftMult
        (
            const scalarField&,
            const SquareMatrix<scalar>&
        );

        const scalarField rightMult
        (
            const SquareMatrix<scalar>&,
            const scalarField&
        );

        SquareMatrix<scalar> outerProd
        (
            const scalarField&,
            const scalarField&
        );

        SquareMatrix<scalar> inv(SquareMatrix<scalar> A);

        //- Compute either global or local sum, based on globalSum flag
        scalar globalSum(const scalarField& field);

        //- Compute either global or local sum, based on globalSum flag
        scalar globalSum(tmp<scalarField>& tfield);

        //- Return optional dictionary with parameters specific to each method
        dictionary coeffsDict();


private:

    // Private Member Functions

        //- No copy construct
        updateMethod(const updateMethod&) = delete;

        //- No copy assignment
        void operator=(const updateMethod&) = delete;


public:

    //- Runtime type information
    TypeName("updateMethod");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            updateMethod,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );


    // Constructors

        //- Construct from components
        updateMethod
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<updateMethod> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~updateMethod() = default;


    // Member Functions

       //- Set objective derivative
       void setObjectiveDeriv(const scalarField& derivs);

       //- Set constraints derivative
       void setConstraintDeriv(const PtrList<scalarField>& derivs);

       //- Set constraints derivative
       void setObjectiveValue(const scalar value);

       //- Set constraints derivative
       void setConstraintValues(const scalarField& values);

       //- Set step for optimisation methods
       void setStep(const scalar eta);

       //- Set globalSum variable.
       //  Should be set by the optimisationType owining the updateMethod
       void setGlobalSum(const bool useGlobalSum);

       //- Return the correction of the design variables
       virtual void computeCorrection()=0;

       //- Return the correction of the design variables
       //const scalarField& returnCorrection() const;

       //- Return the correction of the design variables
       scalarField& returnCorrection();

       void writeCorrection();

       //- Compute merit function. Could be different than the objective
       //- in the presence of constraints
       virtual scalar computeMeritFunction();

       //- Directional derivative of the merit function, in the direction of
       //- the correction. Could be different than the objective directional
       //- derivative in the presence of constraints
       virtual scalar meritFunctionDirectionalDerivative();

       //- Return whether initial eta was set
       bool& initialEtaSet();

       //- Update old correction. useful for quasi-newton methods coupled with
       //- line search
       virtual void updateOldCorrection(const scalarField& oldCorrection);

       //- Write useful quantities to files
       virtual void write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
