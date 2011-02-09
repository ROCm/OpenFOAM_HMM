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
    turbulentTemperatureCoupledMixedSTFvPatchScalarField

Description
    Mixed boundary condition for temperature, to be used for heat-transfer
    on back-to-back baffles.

    Example usage:
        myInterfacePatchName
        {
            type compressible::turbulentTemperatureCoupledMixedST;
            TNbr        T; // name of T field on neighbour region
            K           K; // or none. Name of thermal conductivity field.
                              Calculated from database if equal to 'none'

            QrNbr       Qr; // or none. Name of Qr field on neighbour region
            Qr          Qr; // or none. Name of Qr field on local region

            value       uniform 300;
        }

    Needs to be on underlying directMapped(Wall)FvPatch.

    Note: runs in parallel with arbitrary decomposition. Uses directMapped
    functionality to calculate exchange.

SourceFiles
    turbulentTemperatureCoupledMixedSTFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentTemperatureCoupledMixedSTFvPatchScalarField_H
#define turbulentTemperatureCoupledMixedSTFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "temperatureCoupledBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
    Class turbulentTemperatureCoupledMixedSTFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentTemperatureCoupledMixedSTFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public temperatureCoupledBase
{
    // Private data

        //- Name of field on the neighbour region
        const word TnbrName_;

         //- Name of the radiative heat flux in the neighbout region
        const word QrNbrName_;

        //- Name of the radiative heat flux in local region
        const word QrName_;

        //- Name of thermal conductivity field
        //const word KName_;


public:

    //- Runtime type information
    TypeName("compressible::turbulentTemperatureCoupledMixedST");


    // Constructors

        //- Construct from patch and internal field
        turbulentTemperatureCoupledMixedSTFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentTemperatureCoupledMixedSTFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleMixedFvPatchScalarField onto a
        //  new patch
        turbulentTemperatureCoupledMixedSTFvPatchScalarField
        (
            const
            turbulentTemperatureCoupledMixedSTFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentTemperatureCoupledMixedSTFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        turbulentTemperatureCoupledMixedSTFvPatchScalarField
        (
            const turbulentTemperatureCoupledMixedSTFvPatchScalarField&,
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
                new turbulentTemperatureCoupledMixedSTFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Get corresponding K field
        //tmp<scalarField> K() const;

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
