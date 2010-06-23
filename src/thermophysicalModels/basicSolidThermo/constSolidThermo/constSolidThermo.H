/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    Foam::constSolidThermo

Description
    The thermophysical properties of a constSolidThermo

SourceFiles
    constSolidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef constSolidThermo_H
#define constSolidThermo_H

#include "basicSolidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class constSolidThermo Declaration
\*---------------------------------------------------------------------------*/

class constSolidThermo
:
    public basicSolidThermo
{
        //- Density [kg/m3]
        dimensionedScalar constRho_;

        //- Specific heat capacity [J/(kg.K)]
        dimensionedScalar constCp_;

        //- Thermal conductivity [W/(m.K)]
        //dimensionedSymmTensor constK_;
        dimensionedScalar constK_;

        //- Heat of formation [J/kg]
        dimensionedScalar constHf_;

        //- Emissivity
        dimensionedScalar constEmissivity_;

public:

    //- Runtime type information
    TypeName("constSolidThermo");


    // Constructors

        //- Construct from mesh
        constSolidThermo(const fvMesh& mesh);

    // Destructor

        virtual ~constSolidThermo();


    // Member Functions

        //- Update properties
        virtual void correct();

        //- Density [kg/m3]
        virtual tmp<volScalarField> rho() const;

        //- Specific heat capacity [J/(kg.K)]
        virtual tmp<volScalarField> cp() const;

        //- Thermal conductivity [W/(m.K)]
        //  Note: needs supplied K to be isotropic
        virtual tmp<volScalarField> K() const;

        //- Thermal conductivity [W/(m.K)]
        virtual tmp<volSymmTensorField> directionalK() const;

        //- Heat of formation [J/kg]
        virtual tmp<volScalarField> Hf() const;

        //- Emissivity []
        virtual tmp<volScalarField> emissivity() const;


        // Per patch calculation

            //- Density [kg/m3]
            virtual tmp<scalarField> rho(const label patchI) const;

            //- Specific heat capacity [J/(kg.K)]
            virtual tmp<scalarField> cp(const label patchI) const;

            //- Thermal conductivity [W/(m.K)]
            virtual tmp<scalarField> K(const label patchI) const;

            //- Thermal conductivity [W/(m.K)]
            virtual tmp<symmTensorField> directionalK(const label patchI) const;

            //- Heat of formation [J/kg]
            virtual tmp<scalarField> Hf(const label patchI) const;

            //- Emissivity []
            virtual tmp<scalarField> emissivity(const label) const;


        // I-O

            //- Write the constSolidThermo properties
            virtual bool writeData(Ostream& os) const;

            //- Read solidThermophysicalProperties dictionary
            virtual bool read();

            //- Read solidThermophysicalProperties dictionary
            bool read(const dictionary&);



        // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const constSolidThermo& s);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
