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
    Foam::waveInletVelocityFvPatchVectorField

Description
    Example of the boundary condition specification:
    \verbatim
    inlet
    {
        type            waveInletVelocity;
        value           uniform (0 0 0);
	    leftORright	      1.0;
    }
    \endverbatim

Usage
    \table
        Property    | Description             | Required   | Default value
        type        | type: waveInletVelocity | yes        |
        waveDict    | Dictionary specifying wave variables | no | waveProperties
	    leftORright | Define location of Boundary condition: Left(1) or Right (-1) | yes | -1
    \endtable

Note
    - The value is positive inwards
    - May not work correctly for transonic inlets
    - Strange behaviour with potentialFoam since the U equation is not solved

SourceFiles
    waveInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef waveInletVelocityFvPatchVectorField_H
#define waveInletVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "waveModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
             Class waveInletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class waveInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Dictionary name
        word waveDictName_;

        //- Pointer to the wave model
        tmp<waveModel> waveModel_;

        //- Help function to retrieve the wave model
        waveModel& getWaveModel();


public:

   //- Runtime type information
   TypeName("waveInletVelocity");


   // Constructors

        //- Construct from patch and internal field
        waveInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        waveInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  waveInletVelocityFvPatchVectorField
        //  onto a new patch
        waveInletVelocityFvPatchVectorField
        (
            const waveInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        waveInletVelocityFvPatchVectorField
        (
            const waveInletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new waveInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        waveInletVelocityFvPatchVectorField
        (
            const waveInletVelocityFvPatchVectorField&,
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
                new waveInletVelocityFvPatchVectorField(*this, iF)
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
