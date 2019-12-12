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
    Foam::LBFGS

Description
    The quasi-Newton Limited-memory BFGS formula. Keeps nPrevSteps_ of the
    y and s vectors and approximates the inverse areas through them.
    Values of 3 < nPrevSteps_ < 20 are suggested.

SourceFiles
    LBFGS.C

\*---------------------------------------------------------------------------*/

#ifndef LBFGS_H
#define LBFGS_H

#include "updateMethod.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class LBFGS Declaration
\*---------------------------------------------------------------------------*/

class LBFGS
:
    public updateMethod
{
private:

    // Private Member Functions

        //- No copy construct
        LBFGS(const LBFGS&) = delete;

        //- No copy assignment
        void operator=(const LBFGS&) = delete;


protected:

    // Protected data

        //- Step for the Newton method
        scalar etaHessian_;

        //- Number of first steepest descent steps
        label nSteepestDescent_;

        //- Map to active design variables
        labelList activeDesignVars_;

        //- Number of old corrections and grad differences kept
        label nPrevSteps_;

        //- The previous differences of derivatives. Holds nPrevSteps_ fields
        PtrList<scalarField> y_;

        //- The previous correction. Holds nPrevSteps_ fields
        PtrList<scalarField> s_;

        //- The previous derivatives
        scalarField derivativesOld_;

        //- The previous correction
        scalarField correctionOld_;

        //- Optimisation cycle counter
        label counter_;

        //- Allocate matrices in the first optimisation cycle
        void allocateMatrices();

        //- Move pointers in PtrList to the left and replace last element with
        //- given field
        void pivotFields(PtrList<scalarField>& list, const scalarField& f);

        //- Update y and s vectors
        void updateVectors();

        //- Update design variables
        void update();

        //- Update based on steepest descent
        void steepestDescentUpdate();

        //- Update based on LBFGS
        void LBFGSUpdate();

        //- Read old info from dict
        void readFromDict();


public:

    //- Runtime type information
    TypeName("LBFGS");


    // Constructors

        //- Construct from components
        LBFGS(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~LBFGS() = default;


    // Member Functions

       //- Compute design variables correction
       void computeCorrection();

       //- Update old correction. Useful for quasi-Newton methods coupled with
       //- line search
       virtual void updateOldCorrection(const scalarField& oldCorrection);

       //- Write old info to dict
       virtual void write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
