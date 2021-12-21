/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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
    Foam::compressible::
        turbulentTemperatureRadCoupledMixedFvPatchScalarField

Description
    Mixed boundary condition for temperature and radiation heat transfer
    to be used for in multiregion cases. Optional thin thermal layer
    resistances can be specified through thicknessLayers and kappaLayers
    entries.

    The thermal conductivity \c kappa can either be retrieved from various
    possible sources, as detailed in the class temperatureCoupledBase.

Usage
    \table
        Property     | Description             | Required    | Default value
        Tnbr         | name of the field    | no | T
        qrNbr      | name of the radiative flux in the nbr region | no | none
        qr         | name of the radiative flux in this region | no | none
        thicknessLayers | list of thicknesses per layer [m] | no |
        kappaLayers  | list of thermal conductivites per layer [W/m/K] | no |
        thicknessLayer | single thickness of layer [m] | no |
        kappaLayer  | corresponding thermal conductivity [W/m/K] | no |
        kappaMethod  | inherited from temperatureCoupledBase | inherited |
        kappa        | inherited from temperatureCoupledBase | inherited |
        thermalInertia | Add thermal inertia to wall node | no | false
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            compressible::turbulentTemperatureRadCoupledMixed;
        Tnbr            T;
        qrNbr           qr; // or none. Name of qr field on neighbour region
        qr              qr; // or none. Name of qr field on local region
        thicknessLayers (0.1 0.2 0.3 0.4);
        kappaLayers     (1 2 3 4);
        thermalInertia  false/true;
        kappaMethod     lookup;
        kappa           kappa;
        value           uniform 300;
    }
    \endverbatim

    Needs to be on underlying mapped(Wall)FvPatch.

See also
    Foam::temperatureCoupledBase

SourceFiles
    turbulentTemperatureRadCoupledMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentTemperatureRadCoupledMixedFvPatchScalarField_H
#define turbulentTemperatureRadCoupledMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "temperatureCoupledBase.H"
#include "scalarList.H"
#include "mappedPatchFieldBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
    Class turbulentTemperatureRadCoupledMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentTemperatureRadCoupledMixedFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public temperatureCoupledBase,
    public mappedPatchFieldBase<scalar>
{

    // Private Data

        //- Name of field on the neighbour region
        const word TnbrName_;

        //- Name of the radiative heat flux in the neighbour region
        const word qrNbrName_;

        //- Name of the radiative heat flux in local region
        const word qrName_;

        //- Thickness of layers (either scalarList or a single PatchFunction1)
        scalarList thicknessLayers_;
        autoPtr<PatchFunction1<scalar>> thicknessLayer_;

        //- Conductivity of layers
        scalarList kappaLayers_;
        autoPtr<PatchFunction1<scalar>> kappaLayer_;

        //- Thermal inertia term
        Switch thermalInertia_;


     // Private Functions

            //- Return local alphaSfDelta
            tmp<scalarField> alphaSfDelta() const;

            //- Return delta enthalpy between regions
            tmp<scalarField> deltaH() const;

            //- Return the sum of deltaCoeff*alpha from nbr and local
            tmp<scalarField> beta() const;


            //- Calculate coefficients for assembly matrix
            tmp<Field<scalar>> coeffs
            (
                fvMatrix<scalar>& matrix,
                const Field<scalar>&,
                const label
            ) const;


public:

    //- Runtime type information
    TypeName("compressible::turbulentTemperatureRadCoupledMixed");


    // Constructors

        //- Construct from patch and internal field
        turbulentTemperatureRadCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentTemperatureRadCoupledMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleMixedFvPatchScalarField onto a
        //  new patch
        turbulentTemperatureRadCoupledMixedFvPatchScalarField
        (
            const
            turbulentTemperatureRadCoupledMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentTemperatureRadCoupledMixedFvPatchScalarField
        (
            const turbulentTemperatureRadCoupledMixedFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentTemperatureRadCoupledMixedFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        turbulentTemperatureRadCoupledMixedFvPatchScalarField
        (
            const turbulentTemperatureRadCoupledMixedFvPatchScalarField&,
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
                new turbulentTemperatureRadCoupledMixedFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<scalar>&,
                const labelList&
            );


        //- Given patch temperature calculate corresponding K field. Override
        //- temperatureCoupledBase::kappa to includes effect of any
        //- explicit kappaThickness
        virtual tmp<scalarField> kappa(const scalarField& Tp) const;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Manipulate matrix
        virtual void manipulateMatrix
        (
            fvMatrix<scalar>& m,
            const label iMatrix,
            const direction cmpt
        );


        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
