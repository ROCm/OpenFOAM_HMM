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
    Foam::directionalSolidThermo

Description
    Directional conductivity + table interpolation.

SourceFiles
    directionalSolidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef directionalSolidThermo_H
#define directionalSolidThermo_H

#include "basicSolidThermo.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class directionalSolidThermo Declaration
\*---------------------------------------------------------------------------*/

class directionalSolidThermo
:
    public basicSolidThermo
{
    // Private data

        //- Temperature samples
        Field<scalar> TValues_;

        //- Density at given temperatures
        Field<scalar> rhoValues_;

        Field<scalar> cpValues_;

        Field<vector> KValues_;

        Field<scalar> HfValues_;

        Field<scalar> emissivityValues_;

        //- Coordinate system used for the directional properties
        coordinateSystem coordSys_;

        //- Transformation for cell centres
        volTensorField ccTransforms_;


    // Private Member Functions

        //- Transform principal values of symmTensor
        symmTensor transformPrincipal(const tensor& tt, const vector& st) const;

        //- Transform principal values of symmTensor
        void transformField
        (
            symmTensorField& fld,
            const tensorField& tt,
            const vectorField& st
        ) const;


public:

    //- Runtime type information
    TypeName("directionalSolidThermo");


    // Constructors

        //- Construct from mesh
        directionalSolidThermo(const fvMesh& mesh);

    // Destructor

        virtual ~directionalSolidThermo();


    // Member Functions

        //- Update properties
        virtual void correct();

        //- Density [kg/m3]
        virtual tmp<volScalarField> rho() const;

        //- Specific heat capacity [J/(kg.K)]
        virtual tmp<volScalarField> cp() const;

        //- Thermal conductivity [W/(m.K)]
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
            //  Note: needs Kvalues to be isotropic
            virtual tmp<scalarField> K(const label patchI) const;

            //- Thermal conductivity [W/(m.K)]
            virtual tmp<symmTensorField> directionalK(const label patchI) const;

            //- Heat of formation [J/kg]
            virtual tmp<scalarField> Hf(const label patchI) const;

            //- Emissivity []
            virtual tmp<scalarField> emissivity(const label) const;


        // I-O

            //- Write the directionalSolidThermo properties
            virtual bool writeData(Ostream& os) const;

            //- Read the directionalSolidThermo properties
            virtual bool read();

            //- Read the directionalSolidThermo properties
            bool read(const dictionary& dict);


        // Ostream Operator

        friend Ostream& operator<<
        (
            Ostream& os,
            const directionalSolidThermo& s
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
