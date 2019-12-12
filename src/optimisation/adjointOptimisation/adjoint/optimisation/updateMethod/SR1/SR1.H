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
    Foam::SR1

Description
    The quasi-Newton Symmetric Rank One formula

SourceFiles
    SR1.C

\*---------------------------------------------------------------------------*/

#ifndef SR1_H
#define SR1_H

#include "updateMethod.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                             Class SR1 Declaration
\*---------------------------------------------------------------------------*/

class SR1
:
    public updateMethod
{
protected:

    // Protected data

        //- Step for the Newton method
        scalar etaHessian_;

        //- Number of first steepest descent steps
        label nSteepestDescent_;

        //- For stability check
        scalar ratioThreshold_;

        //- Map to active design variables
        labelList activeDesignVars_;

        //- The Hessian inverse. Should have the size of the active design
        //- variables
        SquareMatrix<scalar> HessianInv_;

        //- The previous Hessian inverse
        SquareMatrix<scalar> HessianInvOld_;

        //- The previous derivatives
        scalarField derivativesOld_;

        //- The previous correction
        scalarField correctionOld_;

        //- Optimisation counter
        label counter_;


    // Protected functions

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
        SR1(const SR1&) = delete;

        //- No copy assignment
        void operator=(const SR1&) = delete;


public:

    //- Runtime type information
    TypeName("SR1");


    // Constructors

        //- Construct from components
        SR1(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~SR1() = default;


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
