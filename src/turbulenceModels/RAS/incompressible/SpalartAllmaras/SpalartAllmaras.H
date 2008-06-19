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
    Foam::incompressible::RAS::SpalartAllmaras

Description
    Spalart-Allmaras 1-eqn mixing-length model for incompressible external
    flows.

    References:
    @verbatim
        "A One-Equation Turbulence Model for Aerodynamic Flows"
        P.R. Spalart,
        S.R. Allmaras,
        La Recherche A´rospatiale, No. 1, 1994, pp. 5–21.

        Extended according to:

        "An Unstructured Grid Generation and Adaptive Solution Technique
        for High Reynolds Number Compressible Flows"
        G.A. Ashford,
        Ph.D. thesis, University of Michigan, 1996.
    @verbatim

SourceFiles
    SpalartAllmaras.C

\*---------------------------------------------------------------------------*/

#ifndef SpalartAllmaras_H
#define SpalartAllmaras_H

#include "RASmodel.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RAS
{

/*---------------------------------------------------------------------------*\
                           Class SpalartAllmaras Declaration
\*---------------------------------------------------------------------------*/

class SpalartAllmaras
:
    public RASmodel
{
    // Private data

        dimensionedScalar alphaNut_;

        dimensionedScalar Cb1_;
        dimensionedScalar Cb2_;
        dimensionedScalar Cw1_;
        dimensionedScalar Cw2_;
        dimensionedScalar Cw3_;
        dimensionedScalar Cv1_;
        dimensionedScalar Cv2_;

        volScalarField nuTilda_;
        volScalarField nut_;

        wallDist d_;


    // Private member functions

        tmp<volScalarField> chi() const;
        tmp<volScalarField> fv1(const volScalarField& chi) const;
        tmp<volScalarField> fv2
        (
            const volScalarField& chi,
            const volScalarField& fv1
        ) const;
        tmp<volScalarField> fv3
        (
            const volScalarField& chi,
            const volScalarField& fv1
        ) const;
        tmp<volScalarField> fw(const volScalarField& Stilda) const;


public:

    //- Runtime type information
    TypeName("SpalartAllmaras");


    // Constructors

        //- from components
        SpalartAllmaras
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport
        );


    // Destructor

        ~SpalartAllmaras()
        {}


    // Member Functions

        //- Return the turbulence viscosity
        tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the effective diffusivity for nuTilda
        tmp<volScalarField> DnuTildaEff() const;

        //- Return the turbulence kinetic energy
        tmp<volScalarField> k() const;

        //- Return the turbulence kinetic energy dissipation rate
        tmp<volScalarField> epsilon() const;

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
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
