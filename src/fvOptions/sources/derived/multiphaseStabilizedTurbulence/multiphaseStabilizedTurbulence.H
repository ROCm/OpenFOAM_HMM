/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd
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
    Foam::fv::multiphaseStabilizedTurbulence

Group
    grpFvOptionsSources

Description
    Applies corrections to the turbulence kinetic energy equation and turbulence
    viscosity field for incompressible multiphase flow cases.

    Turbulence kinetic energy is over-predicted in VOF solvers at the phase
    interface and throughout the water column in nearly-potential flow regions
    beneath surface waves.

    This fvOption applies corrections based on the references:
    \verbatim
        Buoyancy source term in turbulence kinetic energy equation:
            Devolder, B., Rauwoens, P., and Troch, P. (2017).
            Application of a buoyancy-modified k-w SST turbulence model to
            simulate wave run-up around a monopile subjected to regular waves
            using OpenFOAM.
            Coastal Engineering, 125, 81-94.

        Correction to turbulence viscosity:
            Larsen, B.E. and Fuhrman, D.R. (2018).
            On the over-production of turbulence beneath surface waves in
            Reynolds-averaged Navier-Stokes models
            J. Fluid Mech, 853, 419-460
    \endverbatim

Usage
    Example usage:

    \verbatim
    multiphaseStabilizedTurbulence1
    {
        type            multiphaseStabilizedTurbulence;
        active          yes;

        multiphaseStabilizedTurbulenceCoeffs
        {
            // Optional coefficients
            lambda2         0.1;   // A value of 0 sets the nut correction to 0
            Cmu             0.09;  // from k-epsilon model
            C               1.51;  // from k-omega model
            alpha           1.36;  // 1/Prt
        }
    }
    \endverbatim

    The model C coefficient for the k-epsilon model equates to C2/C1 = 1.33;
    the (default) value of 1.51 comes from the k-omega model and is more
    conservative.

SourceFiles
    multiphaseStabilizedTurbulence.C

\*---------------------------------------------------------------------------*/

#ifndef fv_multiphaseStabilizedTurbulence_H
#define fv_multiphaseStabilizedTurbulence_H

#include "fvOption.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class multiphaseStabilizedTurbulence Declaration
\*---------------------------------------------------------------------------*/

class multiphaseStabilizedTurbulence
:
    public option
{
    // Private data

        //- Name of density field
        const word rhoName_;


        // Model coefficients

            //- k-epsilon model Cmu coefficient
            dimensionedScalar Cmu_;

            //- Model coefficient
            //  For k-epsilon models this equates to C2/C1 = 1.33 and for
            //  k-omega = 1.51.  Here the default is the more conservative 1.51
            dimensionedScalar C_;

            //- lambda2 coefficient; default = 0.1
            dimensionedScalar lambda2_;

            //- Buoyancy coefficient
            dimensionedScalar alpha_;


     // Private Member Functions

        //- No copy construct
        multiphaseStabilizedTurbulence
        (
            const multiphaseStabilizedTurbulence&
        ) = delete;

        //- No copy assignment
        void operator=(const multiphaseStabilizedTurbulence&) = delete;


public:

    //- Runtime type information
    TypeName("multiphaseStabilizedTurbulence");


    // Constructors

        //- Construct from explicit source name and mesh
        multiphaseStabilizedTurbulence
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        //- Add explicit contribution to compressible k equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Add explicit contribution to incompressible k equation
        virtual void addSup
        (
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Correct the turbulence viscosity
        virtual void correct(volScalarField& field);

        //- Read source dictionary
        virtual bool read(const dictionary& dict)
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
