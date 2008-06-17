/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::incompressible::RAS::LienCubicKELowRe

Description
    Lien cubic non-linear low-Reynolds k-epsilon turbulence models for
    incompressible flows.

    References:
    @verbatim
        Lien, F.S., Chen, W.L., Leschziner, M.A.,
        "Low-Reynolds-number eddy-viscosity modeling based on non-linear
         stress-strain/vorticity relations"
        Engineering Turbulence Modelling and Experiments 3
        (Edited by Rodi, W. and Bergeles, G.), 91-100. 1996.
        Elsevier Science Publishers.

        Etemad, S., et al.,
        "Turbulent flow and heat transfer in a square-sectioned U bend"
        Progress in compuational fluid dynamics 6, 89-100. 2006.
    @verbatim

SourceFiles
    LienCubicKELowRe.C

\*---------------------------------------------------------------------------*/

#ifndef LienCubicKELowRe_H
#define LienCubicKELowRe_H

#include "turbulenceModel.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RAS
{

/*---------------------------------------------------------------------------*\
                           Class LienCubicKELowRe Declaration
\*---------------------------------------------------------------------------*/

class LienCubicKELowRe
:
    public turbulenceModel
{
    // Private data

        scalar C1;
        scalar C2;
        scalar alphak;
        scalar alphaEps;
        scalar A1;
        scalar A2;
        scalar Ctau1;
        scalar Ctau2;
        scalar Ctau3;
        scalar alphaKsi;

        scalar CmuWall;

        scalar Am;
        scalar Aepsilon;
        scalar Amu;

        volScalarField k_;
        volScalarField epsilon_;

        wallDist y_;

        volTensorField gradU;
        volScalarField eta;
        volScalarField ksi;
        volScalarField Cmu;
        volScalarField fEta;
        volScalarField C5viscosity;

        volScalarField yStar;

        volScalarField nut_;

        volSymmTensorField nonlinearStress;

public:

    //- Runtime type information
    TypeName("LienCubicKELowRe");

    // Constructors

        //- from components
        LienCubicKELowRe
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport
        );


    // Destructor

        ~LienCubicKELowRe()
        {}


    // Member Functions

        //- Return the turbulence viscosity
        tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", alphak*nut_ + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", alphaEps*nut_ + nu())
            );
        }

        //- Return the turbulence kinetic energy
        tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the Reynolds stress tensor
        tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        void correct();

        //- Read turbulenceProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RAS
} // Edn namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
