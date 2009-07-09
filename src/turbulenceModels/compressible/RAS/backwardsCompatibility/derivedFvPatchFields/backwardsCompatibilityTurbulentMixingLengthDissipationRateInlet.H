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
        backwardsCompatibilityTurbulentMixingLengthDissipationRateInletFvPatchScalarField

Description
    Compatibility for the new namespace qualifier
    Foam::compressible::turbulentMixingLengthDissipationRateInletFvPatchScalarField

SourceFiles
    backwardsCompatibilityTurbulentMixingLengthDissipationRateInlet.C

\*---------------------------------------------------------------------------*/

#ifndef backwardsCompatibilityTurbulentMixingLengthDissipationRateInlet_H
#define backwardsCompatibilityTurbulentMixingLengthDissipationRateInlet_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
       Class backwardsCompatibilityTurbulentMixingLengthDissipationRateInletFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class backwardsCompatibilityTurbulentMixingLengthDissipationRateInletFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- turbulent length scale
        scalar mixingLength_;

public:

    //- Runtime type information
    TypeName("turbulentMixingLengthDissipationRateInlet");


    // Constructors

        //- Construct from patch and internal field
        backwardsCompatibilityTurbulentMixingLengthDissipationRateInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        backwardsCompatibilityTurbulentMixingLengthDissipationRateInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  backwardsCompatibilityTurbulentMixingLengthDissipationRateInletFvPatchScalarField
        //  onto a new patch
        backwardsCompatibilityTurbulentMixingLengthDissipationRateInletFvPatchScalarField
        (
            const backwardsCompatibilityTurbulentMixingLengthDissipationRateInletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        backwardsCompatibilityTurbulentMixingLengthDissipationRateInletFvPatchScalarField
        (
            const backwardsCompatibilityTurbulentMixingLengthDissipationRateInletFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new backwardsCompatibilityTurbulentMixingLengthDissipationRateInletFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        backwardsCompatibilityTurbulentMixingLengthDissipationRateInletFvPatchScalarField
        (
            const backwardsCompatibilityTurbulentMixingLengthDissipationRateInletFvPatchScalarField&,
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
                new backwardsCompatibilityTurbulentMixingLengthDissipationRateInletFvPatchScalarField
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
