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
    Foam::RASModels::fWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a wall constraint on the elliptic
    relaxation factor, \c f, which is executed in the \c v2-f based
    turbulence closure models for low- and high-Reynolds number
    turbulent flow cases.

    Reference:
    \verbatim
    Remark on the blending approach for f (tag:PH):
        Popovac, M., & Hanjalić, K. (2007).
        Compound wall treatment for RANS computation of complex
        turbulent flows and heat transfer.
        Flow, turbulence and combustion, 78(2), 177-202.
        DOI:10.1007/s10494-006-9067-x

    Wall-boundary expression for f in kEpsilonPhitF model (tag:LUU):
        Laurence, D. R., Uribe, J. C., & Utyuzhnikov, S. V. (2005).
        A robust formulation of the v2−f model.
        Flow, Turbulence and Combustion, 73(3-4), 169–185.
        DOI:10.1007/s10494-005-1974-8
    \endverbatim

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            fWallFunction;

        // No optional entry
    }
    \endverbatim

Note
    The coefficients \c Cmu, \c kappa, and \c E are obtained from
    the specified \c nutWallFunction in order to ensure that each patch
    possesses the same set of values for these coefficients.

    For \c f, the viscous and inertial sublayer blending approaches were
    claimed to be inviable (PH:p. 194). Therefore, the only blending mode
    for the v2-f model is the stepwise mode where the viscous and inertial
    sublayer contributions switch over a sublayer-intersection value of
    \c y+ estimated from the \c kappa and \c E.

See also
    Foam::fixedValueFvPatchField

SourceFiles
    fWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fWallFunctionFvPatchScalarField_H
#define fWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
               Class fWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchField<scalar>
{
public:

    //- Runtime type information
    TypeName("fWallFunction");


    // Constructors

        //- Construct from patch and internal field
        fWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fWallFunctionFvPatchScalarField
        //- onto a new patch
        fWallFunctionFvPatchScalarField
        (
            const fWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fWallFunctionFvPatchScalarField
        (
            const fWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fWallFunctionFvPatchScalarField
        (
            const fWallFunctionFvPatchScalarField&,
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
                new fWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
