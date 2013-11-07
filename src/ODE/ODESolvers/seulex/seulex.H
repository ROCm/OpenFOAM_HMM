/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::seulex

Description
    An extrapolation-algorithm, based on the linearly implicit Euler method
    with step size control and order selection.

    The implementation is based on the SEULEX code in
    \verbatim
        "Solving Ordinary Differential Equations II: Stiff
         and Differential-Algebraic Problems, second edition",
        Hairer, E.,
        Nørsett, S.,
        Wanner, G.,
        Springer-Verlag, Berlin. 1996.
    \endverbatim

SourceFiles
    seulex.C

\*---------------------------------------------------------------------------*/

#ifndef seulex_H
#define seulex_H

#include "ODESolver.H"

#include "scalarFieldField.H"
#include "scalarMatrices.H"
#include "labelField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class seulex Declaration
\*---------------------------------------------------------------------------*/

class seulex
:
    public ODESolver
{

    // Private data

        //- seulex constants
        static const scalar EPS;

        static const label KMAXX=12, IMAXX = KMAXX + 1;

        static const scalar
            STEPFAC1, STEPFAC2, STEPFAC3, STEPFAC4, STEPFAC5, KFAC1, KFAC2;


        //- Reference to ODESystem
        const ODESystem& ode_;

        // Temporary fields
        mutable label kTarg_;
        mutable labelField nSeq_;
        mutable scalarField cost_, dfdx_, factrl_;
        mutable scalarRectangularMatrix table_, fSave_;
        mutable scalarSquareMatrix dfdy_;
        mutable scalar jacRedo_, theta_;
        mutable bool calcJac_, dense_;
        mutable scalarSquareMatrix a_, coeff_;


        mutable scalar hnext_, hdid_;

        mutable labelList pivotIndices_;

        // Fields space for "step" function
        mutable scalarField hOpt_ ,work_, ySaved_, ySequence_, scale_;

        // Fields used in "dy" function
        mutable scalarField del_, yTemp_, dyTemp_, delTemp_;


    // Private Member Functions

        void step(const scalar& htry, scalar& x,  scalarField& y) const;

        bool dy
        (
            const scalar& x,
            scalarField& y,
            const scalar htot,
            const label k,
            scalarField& yend,
            label& ipt,
            scalarField& scale
        ) const;


        void polyextr
        (
            const label k,
            scalarRectangularMatrix& table,
            scalarField& last
        ) const;


public:

    //- Runtime type information
    TypeName("seulex");


    // Constructors

        //- Construct from ODE
        seulex(const ODESystem& ode, const dictionary& dict);


    // Member Functions
    void solve
    (
        const ODESystem& ode,
        scalar& x,
        scalarField& y,
        scalar& dxTry
    ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
