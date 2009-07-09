/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2006-2009 OpenCFD Ltd.
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
    Foam::compressible::
        backwardsCompatibilityTurbulentMixingLengthFrequencyInletFvPatchScalarField

Description
    Compatibility for the new namespace qualifier
    Foam::compressible::turbulentMixingLengthFrequencyInletFvPatchScalarField

SourceFiles
    backwardsCompatibilityTurbulentMixingLengthFrequencyInletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef compressiblebackwardsCompatibilityTurbulentMixingLengthFrequencyInlet_H
#define compressiblebackwardsCompatibilityTurbulentMixingLengthFrequencyInlet_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
    Class backwardsCompatibilityTurbulentMixingLengthFrequencyInletFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class backwardsCompatibilityTurbulentMixingLengthFrequencyInletFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Turbulent length scale
        scalar mixingLength_;

        //- Name of the turbulent kinetic energy field
        word kName_;


public:

    //- Runtime type information
    TypeName("turbulentMixingLengthFrequencyInlet");


    // Constructors

        //- Construct from patch and internal field
        backwardsCompatibilityTurbulentMixingLengthFrequencyInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        backwardsCompatibilityTurbulentMixingLengthFrequencyInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  backwardsCompatibilityTurbulentMixingLengthFrequencyInletFvPatchScalarField
        //  onto a new patch
        backwardsCompatibilityTurbulentMixingLengthFrequencyInletFvPatchScalarField
        (
            const backwardsCompatibilityTurbulentMixingLengthFrequencyInletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        backwardsCompatibilityTurbulentMixingLengthFrequencyInletFvPatchScalarField
        (
            const backwardsCompatibilityTurbulentMixingLengthFrequencyInletFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new backwardsCompatibilityTurbulentMixingLengthFrequencyInletFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        backwardsCompatibilityTurbulentMixingLengthFrequencyInletFvPatchScalarField
        (
            const backwardsCompatibilityTurbulentMixingLengthFrequencyInletFvPatchScalarField&,
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
                new backwardsCompatibilityTurbulentMixingLengthFrequencyInletFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

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
