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
    Foam::BFGS

Description
    The quasi-Newton BFGS formula

SourceFiles
    BFGS.C

\*---------------------------------------------------------------------------*/

#ifndef BFGS_H
#define BFGS_H

#include "updateMethod.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class BFGS Declaration
\*---------------------------------------------------------------------------*/

class BFGS
:
    public updateMethod
{
private:

    // Private Member Functions

        //- No copy construct
        BFGS(const BFGS&) = delete;

        //- No copy assignment
        void operator=(const BFGS&) = delete;


protected:

    // Protected data

        //- Step for the Newton method
        scalar etaHessian_;

        //- Number of first steepest descent steps
        label nSteepestDescent_;

        //- Map to active design variables
        labelList activeDesignVars_;

        //- Scale the iniitial unitary Hessian approximation
        bool scaleFirstHessian_;

        //- Curvature threshold
        scalar curvatureThreshold_;

        //- The Hessian inverse. Should have the size of the active design
        //- variables
        SquareMatrix<scalar> HessianInv_;

        //- The previous Hessian inverse
        SquareMatrix<scalar> HessianInvOld_;

        //- The previous derivatives
        scalarField derivativesOld_;

        //- The previous correction
        scalarField correctionOld_;

        //- Optimisation cycle counter
        label counter_;

        //- Allocate matrices in the first optimisation cycle
        void allocateMatrices();

        //- Update approximation of the inverse Hessian
        void updateHessian();

        //- Update design variables
        void update();

        //- Read old info from dict
        void readFromDict();


public:

    //- Runtime type information
    TypeName("BFGS");


    // Constructors

        //- Construct from components
        BFGS(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~BFGS() = default;


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
