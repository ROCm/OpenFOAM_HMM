/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016, 2019 OpenFOAM Foundation
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
    Foam::nutWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    The class \c nutWallFunction is a base class that parents the derived
    boundary conditions which provide a wall constraint on various fields, such
    as turbulence kinematic viscosity, i.e. \c nut, for low- and high-Reynolds
    number turbulence models.

    Reference:
    \verbatim
        Default model coefficients:
            Versteeg, H. K., & Malalasekera, W. (2011).
            An introduction to computational fluid dynamics: the finite
            volume method. Harlow: Pearson Education.
            Subsection "3.5.2 k-epsilon model".
    \endverbatim

Usage
    \table
        Property  | Description         | Required   | Default value
        Cmu       | Cmu coefficient     | no         | 0.09
        kappa     | Von Karman constant | no         | 0.41
        E         | E coefficient       | no         | 9.8
    \endtable

    Examples of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            nutWallFunction;

        // Optional entries
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField

SourceFiles
    nutWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutWallFunctionFvPatchScalarField_H
#define nutWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class turbulenceModel;

/*---------------------------------------------------------------------------*\
             Class nutWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
protected:

    // Protected Data

        //- Name of velocity field
        //  Defult is null (not specified) in which case the velocity is
        //  retrieved from the turbulence model
        word UName_;

        //- Cmu coefficient
        scalar Cmu_;

        //- Von Karman constant
        scalar kappa_;

        //- E coefficient
        scalar E_;

        //- Estimated y+ value at the edge of the viscous sublayer
        scalar yPlusLam_;


    // Protected Member Functions

        //- Helper to return the velocity field either from the turbulence
        //- model (default) or the mesh database
        virtual const volVectorField& U(const turbulenceModel& turb) const;

        //- Check the type of the patch
        virtual void checkType();

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const = 0;

        //- Write local wall function variables
        virtual void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("nutWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- nutWallFunctionFvPatchScalarField
        //- onto a new patch
        nutWallFunctionFvPatchScalarField
        (
            const nutWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutWallFunctionFvPatchScalarField
        (
            const nutWallFunctionFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        nutWallFunctionFvPatchScalarField
        (
            const nutWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );


    // Member Functions

        //- Return Cmu
        scalar Cmu() const
        {
            return Cmu_;
        }

        //- Return kappa
        scalar kappa() const
        {
            return kappa_;
        }

        //- Return E
        scalar E() const
        {
            return E_;
        }

        //- Return the nut patchField for the given wall patch
        static const nutWallFunctionFvPatchScalarField& nutw
        (
            const turbulenceModel& turbModel,
            const label patchi
        );

        //- Calculate the y+ at the edge of the viscous sublayer
        static scalar yPlusLam(const scalar kappa, const scalar E);

        //- Return the y+ at the edge of the viscous sublayer
        scalar yPlusLam() const;

        //- Calculate and return the yPlus at the boundary
        virtual tmp<scalarField> yPlus() const = 0;


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
