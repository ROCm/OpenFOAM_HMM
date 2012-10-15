/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::temperatureThermoBaffleFvPatchScalarField

Group
    grpThermoBoundaryConditions

Description
    This boundary condition provides a coupled condition between the primary
    baffle mesh regions.

    The primary region creates it and evolves the thermal baffle heat transfer
    equation.

    The solid thermo and baffle dictionaries are located on the
    primary region.

    type                compressible::temperatureThermoBaffle;

    // Coupled BC.
    neighbourFieldName  T;
    kappa               fluidThermo;
    KName               none;


    // Thermo baffle model
    thermoBaffleModel   thermoBaffle2D;
    regionName          baffleRegion;
    infoOutput          yes;
    active              yes;
    thermoBaffle2DCoeffs
    {
    }


    // Solid thermo
    thermoType
    heSolidThermo
        <pureSolidMixture
            <constIsoSolidTransport
                <constSolidRad
                    <thermo
                        <hConstThermo<incompressible>,sensibleEnthalpy>
                    >
                 >
            >
         >;

    mixture
    {
        specie
        {
            nMoles          1;
            molWeight       20;
        }
        transport
        {
            kappa           0.01;
        }
        radiation
        {
            sigmaS          0;
            kappaRad        0;
            emissivity      1;
        }
        thermodynamics
        {
            Hf              0;
            Cp              15;
        }
        density
        {
            rho             80;
        }
    }

    value               uniform 300;

SeeAlso
    Foam::turbulentTemperatureCoupledBaffleMixedFvPatchScalarField

SourceFiles
    temperatureThermoBaffleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef temperatureThermoBaffleFvPatchScalarField_H
#define temperatureThermoBaffleFvPatchScalarField_H


#include "autoPtr.H"
#include "regionModel.H"
#include "thermoBaffleModel.H"
#include "turbulentTemperatureCoupledBaffleMixedFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
        Class temperatureThermoBaffleFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class temperatureThermoBaffleFvPatchScalarField
:
    public turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
{
    // Private data

        //- Is the baffle owner
        bool owner_;

        //- Thermal baffle
        autoPtr<regionModels::thermoBaffleModels::thermoBaffleModel> baffle_;

        //- Dictionary
        dictionary dict_;


public:

    //- Runtime type information
    TypeName("compressible::temperatureThermoBaffle");


    // Constructors

        //- Construct from patch and internal field
        temperatureThermoBaffleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        temperatureThermoBaffleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  temperatureThermoBaffleFvPatchScalarField onto a new patch
        temperatureThermoBaffleFvPatchScalarField
        (
            const temperatureThermoBaffleFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );


        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new temperatureThermoBaffleFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        temperatureThermoBaffleFvPatchScalarField
        (
            const temperatureThermoBaffleFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new temperatureThermoBaffleFvPatchScalarField(*this, iF)
            );
        }


    // Member functions


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
