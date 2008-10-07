/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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
    Foam::compressible::LESModels::
        muSgsSpalartAllmarasWallFunctionFvPatchScalarField

Description
    Spalart Allmaas wall function boundary condition for compressible flows

SourceFiles
    muSgsSpalartAllmarasWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef muSgsSpalartAllmarasWallFunctionFvPatchScalarField_H
#define muSgsSpalartAllmarasWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
     Class muSgsSpalartAllmarasWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class muSgsSpalartAllmarasWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data


public:

    //- Runtime type information
    TypeName("muSgsSpalartAllmarasWallFunction");


    // Constructors

        //- Construct from patch and internal field
        muSgsSpalartAllmarasWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and Istream
        muSgsSpalartAllmarasWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            Istream&
        );

        //- Construct from patch, internal field and dictionary
        muSgsSpalartAllmarasWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  muSgsSpalartAllmarasWallFunctionFvPatchScalarField
        //  onto a new patch
        muSgsSpalartAllmarasWallFunctionFvPatchScalarField
        (
            const muSgsSpalartAllmarasWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        muSgsSpalartAllmarasWallFunctionFvPatchScalarField
        (
            const muSgsSpalartAllmarasWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new muSgsSpalartAllmarasWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        muSgsSpalartAllmarasWallFunctionFvPatchScalarField
        (
            const muSgsSpalartAllmarasWallFunctionFvPatchScalarField&,
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
                new muSgsSpalartAllmarasWallFunctionFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Evaluation functions

            //- Evaluate the patchField
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::blocking
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
