/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    Foam::compressible::RASModels::mutUWallFunctionFvPatchScalarField

Description
    Wall function boundary condition for walls, based on velocity.

SourceFiles
    mutUWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleMutUWallFunctionFvPatchScalarField_H
#define compressibleMutUWallFunctionFvPatchScalarField_H

#include "mutkWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
            Class mutUWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class mutUWallFunctionFvPatchScalarField
:
    public mutkWallFunctionFvPatchScalarField
{
protected:

    // Protected member functions

        //- Calculate yPLus
        virtual tmp<scalarField> calcYPlus(const scalarField& magUp) const;

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcMut() const;


public:

    //- Runtime type information
    TypeName("mutUWallFunction");


    // Constructors

        //- Construct from patch and internal field
        mutUWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mutUWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  mutUWallFunctionFvPatchScalarField
        //  onto a new patch
        mutUWallFunctionFvPatchScalarField
        (
            const mutUWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mutUWallFunctionFvPatchScalarField
        (
            const mutUWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new mutUWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        mutUWallFunctionFvPatchScalarField
        (
            const mutUWallFunctionFvPatchScalarField&,
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
                new mutUWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Calculate and return the yPlus at the boundary
            virtual tmp<scalarField> yPlus() const;


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
