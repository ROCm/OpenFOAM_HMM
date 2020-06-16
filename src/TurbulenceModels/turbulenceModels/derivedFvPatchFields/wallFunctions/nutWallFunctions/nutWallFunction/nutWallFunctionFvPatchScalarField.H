/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016, 2019 OpenFOAM Foundation
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
    Foam::nutWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    The class \c nutWallFunction is a base class that parents the derived
    boundary conditions which provide a wall constraint on various fields, such
    as turbulent viscosity, i.e. \c nut, or turbulent kinetic energy dissipation
    rate, i.e. \c epsilon, for low- and high-Reynolds number turbulence models.
    The class is not an executable itself, yet a provider for common entries
    to its derived boundary conditions.

    Reference:
    \verbatim
        Default model coefficients (tag:VM):
            Versteeg, H. K., & Malalasekera, W. (2011).
            An introduction to computational fluid dynamics: The finite
            volume method. Harlow: Pearson Education.
            Subsection "3.5.2 k-epsilon model".

        Binomial blending of the viscous and inertial sublayers (tag:ME):
            Menter, F., & Esch, T. (2001).
            Elements of industrial heat transfer prediction.
            In Proceedings of the 16th Brazilian Congress of Mechanical
            Engineering (COBEM), November 2001. vol. 20, p. 117-127.

        Exponential/Max blending of the viscous and inertial sublayers (tag:PH):
            Popovac, M., & Hanjalić, K. (2007).
            Compound wall treatment for RANS computation of complex
            turbulent flows and heat transfer.
            Flow, turbulence and combustion, 78(2), 177-202.
            DOI:10.1007/s10494-006-9067-x
    \endverbatim

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory and other optional entries
        ...

        // Optional (inherited) entries
        Cmu             0.09;
        kappa           0.41;
        E               9.8;
        blending        stepwise;
        n               4.0;
        U               U;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                    | Type   | Req'd | Dflt
      Cmu       | Empirical model coefficient    | scalar | no    | 0.09
      kappa     | von Kármán constant            | scalar | no    | 0.41
      E         | Wall roughness parameter       | scalar | no    | 9.8
      blending  | Viscous/inertial sublayer blending | word | no  | stepwise
      n         | Binomial blending exponent     | scalar | no    | 2.0
      U         | Name of the velocity field     | word   | no    | U
    \endtable

    The inherited entries are elaborated in:
      - \link fixedValueFvPatchFields.H \endlink

    Options for the \c blending entry:
    \verbatim
      stepwise    | Stepwise switch (discontinuous)
      max         | Maximum value switch (discontinuous)
      binomial    | Binomial blending (smooth)
      exponential | Exponential blending (smooth)
    \endverbatim

    wherein \c nut predictions for the viscous and inertial sublayers are
    blended according to the following expressions:

    - \c stepwise (default):

    \f[
        \nu_t = {\nu_t}_{log} \qquad if \quad y^+ > y^+_{lam}
    \f]

    \f[
        \nu_t = {\nu_t}_{vis} \qquad if \quad y^+ <= y^+_{lam}
    \f]

    where
    \vartable
      {\nu_t}_{vis} | \f$\nu_t\f$ prediction in the viscous sublayer
      {\nu_t}_{log} | \f$\nu_t\f$ prediction in the inertial sublayer
      y^+   | estimated wall-normal height of the cell centre in wall units
      y^+_{lam}  | estimated intersection of the viscous and inertial sublayers
    \endvartable


    - \c max (PH:Eq. 27):

    \f[
        \nu_t = max({\nu_t}_{vis}, {\nu_t}_{log})
    \f]


    - \c binomial (ME:Eqs. 15-16):

    \f[
        \nu_t = (({\nu_t}_{vis})^n + ({\nu_t}_{log})^n)^{1/n}
    \f]
    where
    \vartable
      n               | Binomial blending exponent
    \endvartable


    - \c exponential (PH:Eq. 32):

    \f[
        \nu_t = {\nu_t}_{vis} \exp[-\Gamma] + {\nu_t}_{log} \exp[-1/\Gamma]
    \f]

    where (PH:Eq. 31)
    \vartable
      \Gamma       | Blending expression
      \Gamma       | \f$0.01 (y^+)^4 / (1.0 + 5.0 y^+)\f$
    \endvartable

See also
    - Foam::fixedValueFvPatchField

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

    // Protected Enumerations

        //- Options for the blending treatment of viscous and inertial sublayers
        enum blendingType
        {
            STEPWISE,       //!< "Stepwise switch (discontinuous)"
            MAX,            //!< "Maximum value switch (discontinuous)"
            BINOMIAL,       //!< "Binomial blending (smooth)"
            EXPONENTIAL     //!< "Exponential blending (smooth)"
        };

        //- Names for blendingType
        static const Enum<blendingType> blendingTypeNames;


    // Protected Data

        //- Blending treatment (default = blendingType::STEPWISE)
        const enum blendingType blending_;

        //- Binomial blending exponent being used when
        //- blendingType is blendingType::BINOMIAL (default = 4)
        const scalar n_;

        //- Name of velocity field
        //  Default is null (not specified) in which case the velocity is
        //  retrieved from the turbulence model
        word UName_;

        //- Empirical model coefficient
        scalar Cmu_;

        //- von Kármán constant
        scalar kappa_;

        //- Wall roughness parameter
        scalar E_;

        //- Estimated y+ value at the intersection
        //- of the viscous and inertial sublayers
        scalar yPlusLam_;


    // Protected Member Functions

        //- Helper to return the velocity field either from the turbulence
        //- model (default) or the mesh database
        virtual const volVectorField& U(const turbulenceModel& turb) const;

        //- Check the type of the patch
        virtual void checkType();

        //- Calculate the turbulent viscosity
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

        // No clone methods - abstract class


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

        //- Estimate the y+ at the intersection of the two sublayers
        static scalar yPlusLam(const scalar kappa, const scalar E);

        //- Return the estimated y+ at the two-sublayer intersection
        scalar yPlusLam() const;

        //- Return the blended nut according to the chosen blending treatment
        scalar blend
        (
            const scalar nutVis,
            const scalar nutLog,
            const scalar yPlus
        ) const;

        //- Calculate and return the yPlus at the boundary
        //  yPlus is the first-cell-centre height from boundary in wall units
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
