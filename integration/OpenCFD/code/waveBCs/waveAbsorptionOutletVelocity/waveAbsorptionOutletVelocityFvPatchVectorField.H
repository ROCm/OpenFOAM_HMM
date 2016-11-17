/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::waveAbsorptionOutletVelocityFvPatchVectorField

Description
    Example of the boundary condition specification:
    \verbatim
    inlet
    {
        type            waveAbsorptionOutletVelocity;
        waveDict        wavesDict;
        value           uniform (0 0 0);
    }
    \endverbatim

Usage
    \table
        Property     | Description             | Required    | Default value
        type         | type: waveAbsorptionOutletVelocity | yes         |
        waveDict    | Dictionary specifying wave variables | no | waveProperties
    \endtable

Note
    - The value is positive inwards
    - May not work correctly for transonic inlets
    - Strange behaviour with potentialFoam since the U equation is not solved

SourceFiles
    waveAbsorptionOutletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef waveAbsorptionOutletVelocityFvPatchVectorField_H
#define waveAbsorptionOutletVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
       Class waveAbsorptionOutletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class waveAbsorptionOutletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Dictionary name
        word waveDictName_;


public:

   //- Runtime type information
   TypeName("waveAbsorptionOutletVelocity");


   // Constructors

        //- Construct from patch and internal field
        waveAbsorptionOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        waveAbsorptionOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  waveAbsorptionOutletVelocityFvPatchVectorField
        //  onto a new patch
        waveAbsorptionOutletVelocityFvPatchVectorField
        (
            const waveAbsorptionOutletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        waveAbsorptionOutletVelocityFvPatchVectorField
        (
            const waveAbsorptionOutletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new waveAbsorptionOutletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        waveAbsorptionOutletVelocityFvPatchVectorField
        (
            const waveAbsorptionOutletVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new waveAbsorptionOutletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

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
