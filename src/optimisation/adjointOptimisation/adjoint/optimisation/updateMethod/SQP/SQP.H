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
    Foam::SQP

Description
    The quasi-Newton SQP formula for constrained optimisation

SourceFiles
    SQP.C

\*---------------------------------------------------------------------------*/

#ifndef SQP_H
#define SQP_H

#include "constrainedOptimisationMethod.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class SQP Declaration
\*---------------------------------------------------------------------------*/

class SQP
:
    public constrainedOptimisationMethod
{
protected:

    // Protected data

        //- Step for the Newton method
        scalar etaHessian_;

        //- Map to active design variables
        labelList activeDesignVars_;

        //- Scale the initial unitary Hessian approximation
        bool scaleFirstHessian_;

        //- Curvature threshold
        scalar dumpingThreshold_;

        //- Derivatives of the Lagrangian function
        scalarField LagrangianDerivatives_;

        //- The Hessian inverse. Should have the size of the active design
        //- variables
        SquareMatrix<scalar> Hessian_;

        //- The previous Hessian inverse
        SquareMatrix<scalar> HessianOld_;

        //- The previous objective derivatives
        scalarField objectiveDerivativesOld_;

        //- The previous constraint derivatives
        List<scalarField> constraintDerivativesOld_;

        //- The previous correction
        scalarField correctionOld_;

        //- Lagrange multipliers
        scalarField lamdas_;

        //- Optimisation cycle count
        label counter_;

        //- Name of the objective folder
        fileName objFunctionFolder_;

        //- File including the l1 merit function
        autoPtr<OFstream> meritFunctionFile_;

        //- Penalty value for the merit function
        scalar mu_;

        //- Safety factor
        scalar delta_;


private:

    // Private Member Functions

        //- No copy construct
        SQP(const SQP&) = delete;

        //- No copy assignment
        void operator=(const SQP&) = delete;

        //- Make folder holding the Lagrangian file
        void makeFolder();

        //- Allocate fields and matrices when size of design variables
        //- is known
        void allocateMatrices();

        //- Update the Hessian matrix
        void updateHessian();

        //- Compute new lamdas and update correction
        void computeLagrangeMultipliersAndCorrect();

        //- Store old fields needed for the next iter
        void storeOldFields();

        //- Read old values from dict, if present
        void readFromDict();


public:

    //- Runtime type information
    TypeName("SQP");


    // Constructors

        //- Construct from components
        SQP(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~SQP() = default;


    // Member Functions

       //- Compute design variables correction
       void computeCorrection();

       //- Compute merit function. Could be different than the objective
       //- in the presence of constraints
       virtual scalar computeMeritFunction();

       //- Derivative of the merit function. Could be different than the
       //- objective derivative in the presence of constraints
       virtual scalar meritFunctionDirectionalDerivative();

       //- Update old correction. Useful for quasi-Newton methods coupled with
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
