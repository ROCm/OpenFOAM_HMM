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

        // SEULEX constants
        static const label kMaxx_ = 12;
        static const label iMaxx_ = kMaxx_ + 1;

        static const scalar
            stepFactor1_, stepFactor2_, stepFactor3_,
            stepFactor4_, stepFactor5_,
            kFactor1_, kFactor2_;

        // Temporary fields
        mutable label kTarg_;
        mutable labelField nSeq_;
        mutable scalarField cost_, dfdx_, factrl_;
        mutable scalarRectangularMatrix table_, fSave_;
        mutable scalarSquareMatrix dfdy_;
        mutable scalar jacRedo_, theta_;
        mutable bool calcJac_, dense_;
        mutable scalarSquareMatrix a_, coeff_;
        mutable labelList pivotIndices_;

        // Fields space for "solve" function
        mutable scalarField dxOpt_, work_, y0_, ySequence_, scale_;

        // Fields used in "dy" function
        mutable scalarField del_, yTemp_, dyTemp_, delTemp_;


    // Private Member Functions

        bool dy
        (
            const scalar x,
            scalarField& y,
            const scalar dxTot,
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
