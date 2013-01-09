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
    Foam::thermoBaffleFvPatchScalarField

Group
    grpThermoBoundaryConditions

Description
    This boundary condition provides a coupled temperature condition between
    multiple mesh regions.  The regions are generally referred to as the:
    - primary region, and
    - baffle region

    The primary region creates the baffle region and evolves its energy
    equation either:
    - 1-D, normal to each patch face
    - 2-D, normal and tangential components

    The thermodynamic properties of the baffle material are specified via
    dictionary entries.

    \heading Patch usage

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type                compressible::temperatureThermoBaffle;

        // Coupled boundary condition
        neighbourFieldName  T;
        kappa               fluidThermo;
        KName               none;


        // Thermo baffle model
        thermoBaffleModel   thermoBaffle;
        regionName          baffleRegion;
        infoOutput          yes;
        active              yes;
        thermoBaffleCoeffs
        {
        }


        // Solid thermo
        thermoType
        {
            type            heSolidThermo;
            mixture         pureSolidMixture;
            transport       constIso;
            thermo          hConst;
            equationOfState rhoConst;
            specie          specie;
            energy          sensibleEnthalpy;
        }

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

        radiation
        {
            radiationModel  opaqueSolid;
            absorptionEmissionModel none;
            scatterModel    none;
        }

        value               uniform 300;
    }
    \endverbatim

SeeAlso
    Foam::turbulentTemperatureCoupledBaffleMixedFvPatchScalarField
    Foam::regionModels::thermoBaffleModels::thermoBaffleModel

SourceFiles
    thermoBaffleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef thermoBaffleFvPatchScalarField_H
#define thermoBaffleFvPatchScalarField_H


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
        Class thermoBaffleFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class thermoBaffleFvPatchScalarField
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
    TypeName("compressible::thermoBaffle");


    // Constructors

        //- Construct from patch and internal field
        thermoBaffleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        thermoBaffleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  thermoBaffleFvPatchScalarField onto a new patch
        thermoBaffleFvPatchScalarField
        (
            const thermoBaffleFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );


        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new thermoBaffleFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        thermoBaffleFvPatchScalarField
        (
            const thermoBaffleFvPatchScalarField&,
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
                new thermoBaffleFvPatchScalarField(*this, iF)
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
