/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenCFD Ltd
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

    Foam::
    compressible::
    humidityTemperatureCoupledMixedFvPatchScalarField

Description
    Mixed boundary condition for temperature to be used on coupling flow and
    solid regions. This BC can operate in four modes:

    1) 'inert' : thermal inertia is important and no condensation/evaporation
                  is taken place.
    2) 'condensation' : just condensation is taken place
    3) 'vaporization' : just evaporation is take place
    4) 'condEvap' : both condensation and evaporation take place

    For 'inert' operation the 'rho', 'thickness' and 'cp' entries are needed.

    In 'condensation' mode when the wall temperature (Tw) is bellow the dew
    temperature (Tdew) condesation takes place and the resulting condensed mass
    is stored on the wall.

    In 'vaporization' the initial mass is vaporized when Tw is above the
    input vaporization temperature (Tvap).

    In 'condEvap', condensation and evaporation take place simultaneously.

    The BC assumes no mass flow on the wall.i.e the mass condensed on a face
    remains on that face. It uses a 'lump mass' model to include thermal
    inertia effects.

    It assumes a drop-wise type of condensation and its heat transfer Nu number
    is:

        51104 + 2044*T   T > 295  T < 373
        255510           T > 373

    T. Bergam, A.Lavine, F. Incropera and D. Dewitt. Heat and Mass Transfer.
    7th Edition. Chapter 10.

    The mass transfer correlation used is hm = Dab*Sh/L

    where:

        Dab is the mass vapor difussivity
        L   is the characteristic lenght
        Sc  the Schmidt number and it is calculated as:

            0.664*sqrt(Re)*cbrt(Sc)  Re < 5.0E+05
            0.037*pow(Re, 0.8)*cbrt(Sc)  Re > 5.0E+05

    NOTE: The correclation used to calculate Tdew is for water vapor.
    In addition a scalar transport for the carrier specie have to be specified
    via function objects or in the main solver. This specie transports the
    vapour phase in the main ragion. The BC of this specie on the coupled wall
    has to fixedGradient in order to allow condensation or evaporation of the
    vapor in or out of this wall


    Example usage:

    On the fluid side

        myInterfacePatchName
        {
            type            thermalHumidityCoupledMixed;
            kappa           fluidThermo;
            kappaName       none;

            // Modes of operation: inert, condensation, vaporization, condEvap
            mode            condEvap;

            // Carrier species name
            specieName      H2O;

            // Carrier molecular weight
            carrierMolWeight           28.9;

            // Characteristic lenght of the wall
            L               0.1;

            // Vaporasation temperature
            Tvap            273;

            // Liquid properties for the condensed mass
            liquid
            {
                H2O
                {
                    defaultCoeffs       yes;
                }
            }

            // thickness, density and cp required for inert and condensation
            // modes

            //thickness       uniform 0;
            //cp              uniform 0;
            //rho             uniform 0;

            value           $internalField;
        }


    On the solid side:

        myInterfacePatchName
        {
            type            thermalInertiaMassTransferCoupledMixed;
            kappa           solidThermo;
            kappaName       none;
            value           uniform 260;
        }


SourceFiles
    humidityTemperatureCoupledMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef humidityTemperatureCoupledMixedFvPatchScalarField_H
#define humidityTemperatureCoupledMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "temperatureCoupledBase.H"
#include "liquidProperties.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class humidityTemperatureCoupledMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class humidityTemperatureCoupledMixedFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public temperatureCoupledBase
{
public:

    //  Public enumeration

        //- Modes of mass transfer
        enum massTransferMode
        {
            mConstantMass,
            mCondensation,
            mEvaporation,
            mCondensationAndEvaporation
        };


private:

    // Private data

        static const NamedEnum<massTransferMode, 4> MassModeTypeNames_;

        //- BC mode
        massTransferMode mode_;

        //- Name of field on the neighbour region
        const word TnbrName_;

         //- Name of the radiative heat flux in the neighbout region
        const word QrNbrName_;

        //- Name of the radiative heat flux
        const word QrName_;

        //- Name of the species on which the mass transfered (default H2O)
        const word specieName_;

        //- Liquid properties
        autoPtr<liquidProperties> liquid_;

        //- Liquid dictionary
        dictionary liquidDict_;

        //- Mass accumulated on faces
        scalarField mass_;

        //- Vaporization temperature
        scalar Tvap_;

        //- Cache myDelta
        scalarField myKDelta_;

        //- Phase change energy
        scalarField dmHfg_;

        //- Thermal inertia
        scalarField mpCpTp_;

        //- Average molecular weight for the carrier mixture in the gas phase
        scalar Mcomp_;

        //- Characteristic length scale
        scalar L_;

        //- Fluid side
        bool fluid_;

        //- Cp field for inert mode
        scalarField cp_;

        //- Thickness field for inert mode
        scalarField thickness_;

        //- Density field for inert mode
        scalarField rho_;


    // Private members

        //- Calculation of Sh
        scalar Sh(const scalar Re,const scalar Sc) const;

        //- Calculation of htc from the condensed surface
        scalar htcCondensation(const scalar TSat,const scalar Re) const;

        //- Create thickness field for output
        volScalarField& thicknessField(const word& , const fvMesh&);


public:

    //- Runtime type information
    TypeName("humidityTemperatureCoupledMixed");


    // Constructors

        //- Construct from patch and internal field
        humidityTemperatureCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        humidityTemperatureCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleMixedFvPatchScalarField onto a
        //  new patch
        humidityTemperatureCoupledMixedFvPatchScalarField
        (
            const
            humidityTemperatureCoupledMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new humidityTemperatureCoupledMixedFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        humidityTemperatureCoupledMixedFvPatchScalarField
        (
            const humidityTemperatureCoupledMixedFvPatchScalarField&,
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
                new humidityTemperatureCoupledMixedFvPatchScalarField
                (
                    *this,
                    iF
                )
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


        //- Return myKDelta
        const scalarField myKDelta() const
        {
            return myKDelta_;
        }

        //- Return mpCpdTpd
        const scalarField mpCpTp() const
        {
            return mpCpTp_;
        }

        //- Return dmHfg
        const scalarField dmHfg() const
        {
            return dmHfg_;
        }


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
