/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::variableHeightFlowRateFvPatchScalarField

Description
    This boundary condition uses zeroGradient within a specified range of
    values for phase fraction alpha. The range is defined within the
    boundary condition by the lowerBound and upperBound.

    alpha > upperBound: fixedValue with uniform value of upperBound
    lowerBound <= alpha <= upperBound: zeroGradient
    alpha < lowerBound: fixedValue with uniform value of lowerBound

    Example:
    \verbatim
        inlet
        {
            type            clippedZeroGradient;
            lowerBound      0.0;
            upperBound      0.9;
            value           uniform 0;
        }
    \verbatim

SourceFiles
    variableHeightFlowRateFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef variableHeightFlowRateFvPatchScalarField_H
#define variableHeightFlowRateFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class variableHeightFlowRateFvPatchScalar Declaration
\*---------------------------------------------------------------------------*/

class variableHeightFlowRateFvPatchScalarField
:
    public mixedFvPatchScalarField
{

protected:

    // Protected data

        //- Name of flux field
        word phiName_;

        //- Lower bound for alpha1
        scalar lowerBound_;

        //- Upper bound for alpha1
        scalar upperBound_;

public:

    //- Runtime scalar information
    TypeName("variableHeightFlowRate");


    // Constructors

        //- Construct from patch and internal field
        variableHeightFlowRateFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        variableHeightFlowRateFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  variableHeightFlowRateFvPatchScalarField onto a new patch
        variableHeightFlowRateFvPatchScalarField
        (
            const variableHeightFlowRateFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        variableHeightFlowRateFvPatchScalarField
        (
            const variableHeightFlowRateFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar> > clone() const
        {
            return tmp<fvPatchField<scalar> >
            (
                new variableHeightFlowRateFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        variableHeightFlowRateFvPatchScalarField
        (
            const variableHeightFlowRateFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar> > clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar> >
            (
                new variableHeightFlowRateFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
