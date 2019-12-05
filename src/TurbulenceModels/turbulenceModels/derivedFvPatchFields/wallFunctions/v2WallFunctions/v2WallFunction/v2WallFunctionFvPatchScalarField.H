/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016, 2019 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::RASModels::v2WallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a wall constraint on wall-normal velocity
    scale, i.e. \c v2, for low- and high-Reynolds number turbulence models.

    The model operates in two modes, based on the computed viscous-to-turbulent
    switch-over \c y+ value derived from \c kappa and \c E.

Usage
    \table
        Property     | Description             | Required    | Default value
        Cv2          | model coefficient       | no          |  0.193
        Bv2          | model coefficient       | no          | -0.94
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            v2WallFunction;

        // Optional entries
    }
    \endverbatim

Note
    The coefficients \c Cmu, \c kappa, and \c E are obtained from
    the specified \c nutWallFunction in order to ensure that each patch
    possesses the same set of values for these coefficients.

See also
    Foam::fixedValueFvPatchField

SourceFiles
    v2WallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef v2WallFunctionFvPatchScalarField_H
#define v2WallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
              Class v2WallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class v2WallFunctionFvPatchScalarField
:
    public fixedValueFvPatchField<scalar>
{
protected:

    // Protected Data

        //- Cv2 coefficient
        scalar Cv2_;

        //- Bv2 coefficient
        scalar Bv2_;


public:

    //- Runtime type information
    TypeName("v2WallFunction");


    // Constructors

        //- Construct from patch and internal field
        v2WallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        v2WallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given v2WallFunctionFvPatchScalarField
        //- onto a new patch
        v2WallFunctionFvPatchScalarField
        (
            const v2WallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        v2WallFunctionFvPatchScalarField
        (
            const v2WallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new v2WallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        v2WallFunctionFvPatchScalarField
        (
            const v2WallFunctionFvPatchScalarField&,
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
                new v2WallFunctionFvPatchScalarField(*this, iF)
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

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
