/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016, 2019 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::kLowReWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a wall constraint on the turbulent kinetic
    energy, i.e. \c k, for low- and high-Reynolds number turbulence models.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries (unmodifiable)
        type            kLowReWallFunction;

        // Optional entries (unmodifiable)
        Ceps2           1.9;
        Ck              -0.416;
        Bk              8.366;
        C               11.0;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                     | Type   | Req'd | Dflt
      type      | Type name: kLowReWallFunction   | word   | yes   | -
      Ceps2     | Model coefficient               | scalar | no    |  1.9
      Ck        | Model coefficient               | scalar | no    | -0.416
      Bk        | Model coefficient               | scalar | no    |  8.366
      C         | Model coefficient               | scalar | no    |  11.0
    \endtable

    The inherited entries are elaborated in:
      - \link fixedValueFvPatchField.H \endlink

    Viscous and inertial sublayer predictions for \c k are blended in
    a stepwise manner:

    \f[
        k = k_{log} \qquad if \quad y^+ > y^+_{lam}
    \f]
    \f[
        k = k_{vis} \qquad if \quad y^+ <= y^+_{lam}
    \f]
    where
    \vartable
      k_{vis}   | k prediction in the viscous sublayer
      k_{log}   | k prediction in the inertial sublayer
      y^+       | estimated wall-normal height of the cell centre in wall units
      y^+_{lam} | estimated intersection of the viscous and inertial sublayers
    \endvartable

Note
    The coefficients \c Cmu, \c kappa, and \c E are obtained from
    the specified \c nutWallFunction in order to ensure that each patch
    possesses the same set of values for these coefficients.

See also
    - Foam::fixedValueFvPatchField
    - Foam::kqRWallFunctionFvPatchScalarField

SourceFiles
    kLowReWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef kLowReWallFunctionFvPatchScalarField_H
#define kLowReWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class kLowReWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class kLowReWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchField<scalar>
{
protected:

    // Protected Data

        //- Ceps2 coefficient
        scalar Ceps2_;

        //- Ck coefficient
        scalar Ck_;

        //- Bk coefficient
        scalar Bk_;

        //- C coefficient
        scalar C_;


public:

    //- Runtime type information
    TypeName("kLowReWallFunction");


    // Constructors

        //- Construct from patch and internal field
        kLowReWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        kLowReWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given kLowReWallFunctionFvPatchScalarField
        //- onto a new patch
        kLowReWallFunctionFvPatchScalarField
        (
            const kLowReWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        kLowReWallFunctionFvPatchScalarField
        (
            const kLowReWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new kLowReWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        kLowReWallFunctionFvPatchScalarField
        (
            const kLowReWallFunctionFvPatchScalarField&,
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
                new kLowReWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
