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
    Foam::DBFGS

Description
    The quasi-Newton BFGS formula with the dampening proposed by Powell

SourceFiles
    DBFGS.C

\*---------------------------------------------------------------------------*/

#ifndef DBFGS_H
#define DBFGS_H

#include "updateMethod.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class DBFGS Declaration
\*---------------------------------------------------------------------------*/

class DBFGS
:
    public updateMethod
{
protected:

    // Protected data

        //- Step for the Newton method
        scalar etaHessian_;

        //- Number of first steepest descent steps
        label nSteepestDescent_;

        //- Map to active design variables
        labelList activeDesignVars_;

        //- Scale the initial unitary Hessian approximation
        bool scaleFirstHessian_;

        //- Curvature threshold
        scalar curvatureThreshold_;

        //- The Hessian. Should have the size of the active design variables
        SquareMatrix<scalar> Hessian_;

        //- The previous Hessian
        SquareMatrix<scalar> HessianOld_;

        //- The previous derivatives
        scalarField derivativesOld_;

        //- The previous correction
        scalarField correctionOld_;

        //- Optimisation cycle counter
        label counter_;

        //- Threshold for damping
        scalar gamma_;


    // Protected Member Functions

        //- Allocate matrices in the first optimisation cycle
        void allocateMatrices();

        //- Update approximation of the inverse Hessian
        void updateHessian();

        //- Update design variables
        void update();

        //- Read old info from dict
        void readFromDict();

private:

    // Private Member Functions

        //- No copy construct
        DBFGS(const DBFGS&) = delete;

        //- No copy assignment
        void operator=(const DBFGS&) = delete;


public:

    //- Runtime type information
    TypeName("DBFGS");

    // Constructors

        //- Construct from components
        DBFGS(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~DBFGS() = default;


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
