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
    Foam::RKF45

Description
    4/5th Order Runge-Kutta-Fehlberg ODE solver

    References:
    \verbatim
        "Low-order classical Runge-Kutta formulas with step size control
         and their application to some heat transfer problems."
         Fehlberg, E.,
         NASA Technical Report 315 1969

         "Solving Ordinary Differential Equations I: Nonstiff Problems,
          second edition",
         Hairer, E.,
         Nørsett, S.,
         Wanner, G.,
         Springer-Verlag, Berlin. 1993, ISBN 3-540-56670-8.
    \endverbatim

    This method embedds the 4-th order integration step into the 5-th order step
    and allows to perform an adapdive step-size control using these two order
    without the need of re-evaluation.

    \f{align}{
     k_1 &= h f(t_k, y_k) \\
     k_2 &= h f(t_k + \frac{1}{4}h, y_k + \frac{1}{4}k_1) \\
     k_3 &= h f(t_k + \frac{3}{8}h, y_k + \frac{3}{32}k_1 + \frac{9}{32}k_2) \\
     k_4 &= h f(t_k + \frac{12}{13}h, y_k + \frac{1932}{2197}k_1
              \frac{7200}{2197}k_2 + \frac{7296}{2197}k_3) - \\
     k_5 &= h f(t_k + h, y_k + \frac{439}{216}k_1 - 8k_2 + \frac{3680}{513}k_3 -
           \frac{845}{4104}k_4) \\
     k_6 &= h f(t_k + \frac{1}{2}h, y_k - \frac{8}{27}k_1 + 2k_2 -
           \frac{3544}{2565}k_3 + \frac{1859}{4104}k_4 - \frac{11}{40}k_5) \\
    \f}

    Which yields the following update-steps for the 4th and 5th order:

    \f{align}{
     \Delta_4 &= \frac{25}{216}k_1 + \frac{1408}{2565}k_3 +
                 \frac{2197}{4101}k_4 - \frac{1}{5}k_5 \\
     \Delta_5 &= \frac{16}{135}k_1 + \frac{6656}{12825}k_3 +
                 \frac{9}{50}k_5 + \frac{2}{55}k_6
    \f}

    The difference between the 5th and 4th order (\f$\epsilon = \left|\Delta_5
    - \Delta_4\right|\f$) can be used to determine if the stepsize \f$h\f$ needs
    to be adjusted. The step-size does not need to be reduced if the following
    inequation holds true:

    \f[
      \left|\Delta_5 - \Delta_4\right| < \epsilon_{abs} +
                                         \epsilon_{rel}\left|y\right|
    \f]

    where \f$\left|y\right|\f$ is the element wise maximum of
    \f$\left|y_k\right|\f$ and \f$\left|y_{k+1}\right|\f$.

SourceFiles
    RKF45.C

\*---------------------------------------------------------------------------*/

#ifndef RKF45_H
#define RKF45_H

#include "ODESolver.H"
#include "adaptiveSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class RKF45 Declaration
\*---------------------------------------------------------------------------*/

class RKF45
:
    public ODESolver,
    public adaptiveSolver
{
    // Private data

        //- RKF45 Constants
        static const scalar
            c2, c3, c4, c5, c6,
            a21, a31, a32, a41, a42, a43, a51, a52, a53, a54,
            a61, a62, a63, a64, a65,
            b41, b42, b43, b44, b45, b46,
            b51, b52, b53, b54, b55, b56;


        // Temporary fields
        mutable scalarField yTemp_;
        mutable scalarField k2_;
        mutable scalarField k3_;
        mutable scalarField k4_;
        mutable scalarField k5_;
        mutable scalarField k6_;

        //- Error-estimate field
        mutable scalarField err_;


public:

    //- Runtime type information
    TypeName("RKF45");


    // Constructors

        //- Construct from ODE
        RKF45(const ODESystem& ode, const dictionary& dict);


    // Member Functions

        //- Solve a single step dx and return the error
        scalar solve
        (
            const ODESystem& ode,
            const scalar x0,
            const scalarField& y0,
            const scalarField& dydx0,
            const scalar dx,
            scalarField& y
        ) const;

        //- Solve the ODE system and the update the state
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
