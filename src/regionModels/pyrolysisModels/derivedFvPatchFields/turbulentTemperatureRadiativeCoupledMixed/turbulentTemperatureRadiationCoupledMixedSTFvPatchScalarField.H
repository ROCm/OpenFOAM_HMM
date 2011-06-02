/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2011 OpenCFD Ltd.
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
    turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField

Description
    Mixed boundary condition for temperature, to be used for heat-transfer
    on back-to-back baffles.

    If my temperature is T1, neighbour is T2:

    T1 > T2: my side becomes fixedValue T2 bc, other side becomes fixedGradient.


    Example usage:
        myInterfacePatchName
        {
            type compressible::turbulentTemperatureRadiationCoupledMixedST;
            neighbourFieldName  T;
            K                   lookup;
            KName               K;

            neighbourFieldRadiativeName Qr; // or none. Use on solid side
            fieldRadiativeName  Qr;         // fluid. Use on fluid side

            value       uniform 300;
        }

    Needs to be on underlying directMapped(Wall)FvPatch.

    Note: K : heat conduction at patch. Gets supplied how to lookup/calculate K:
    - 'lookup' : lookup volScalarField (or volSymmTensorField) with name
    - 'basicThermo' : use basicThermo and compressible::RASmodel to calculate K
    - 'solidThermo' : use basicSolidThermo K()
    - 'directionalSolidThermo' directionalK()

    Note: runs in parallel with arbitrary decomposition. Uses directMapped
    functionality to calculate exchange.

    Note: lags interface data so both sides use same data.
    - problem: schedule to calculate average would interfere
    with standard processor swaps.
    - so: updateCoeffs sets both to same Twall. Only need to do
    this for last outer iteration but don't have access to this.

SourceFiles
    turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField_H
#define turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField_H

#include "mixedFvPatchFields.H"
//#include "temperatureCoupledBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
   Class turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField
   Declaration
\*---------------------------------------------------------------------------*/

class turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField
:
    public mixedFvPatchScalarField
//    public temperatureCoupledBase
{
    // Private data

        //- Name of field on the neighbour region
        const word neighbourFieldName_;

         //- Name of the radiative heat flux in the neighbout region
        const word neighbourFieldRadiativeName_;

        //- Name of the radiative heat flux in the my region
        const word fieldRadiativeName_;

        //- Name of thermal conductivity field
        const word KName_;

        //- how to obtain radiative flux
        enum operationMode
        {
            radFromNeighbour,
            radFromMe,
            noRad,
            unknown
        };
        static const NamedEnum<operationMode, 4> operationModeNames;

        //- Previous iteration mode
        operationMode oldMode_;


public:

    //- Runtime type information
    TypeName("compressible::turbulentTemperatureRadiationCoupledMixedST");


    // Constructors

        //- Construct from patch and internal field
        turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleMixedFvPatchScalarField onto a
        //  new patch
        turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField
        (
            const
            turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new
                turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField
        (
            const 
                turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField&,
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
                new 
                turbulentTemperatureRadiationCoupledMixedSTFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Get corresponding K field
        tmp<scalarField> K() const;

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
