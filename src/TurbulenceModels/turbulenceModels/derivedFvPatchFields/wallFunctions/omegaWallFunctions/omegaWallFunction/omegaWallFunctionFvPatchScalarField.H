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
    Foam::omegaWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a wall constraint on the specific
    dissipation rate, i.e. \c omega, and the turbulent kinetic energy
    production contribution, i.e. \c G, for low- and high-Reynolds number
    turbulence models.

    Reference:
    \verbatim
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

        Tanh blending of the viscous and inertial sublayers (tag:KAS):
            Knopp, T., Alrutz, T., & Schwamborn, D. (2006).
            A grid and flow adaptive wall-function method for RANS
            turbulence modelling.
            Journal of Computational Physics, 220(1), 19-40.
            DOI:10.1016/j.jcp.2006.05.003
    \endverbatim

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries (unmodifiable)
        type            omegaWallFunction;

        // Optional entries (unmodifiable)
        beta1           0.075;
        blending        binomial2;
        n               2.0;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    \table
      Property  | Description                     | Type   | Req'd  | Dflt
      type      | Type name: omegaWallFunction    | word   | yes    | -
      beta1     | Model coefficient               | scalar | no     | 0.075
      blending  | Viscous/inertial sublayer blending method <!--
                                              --> | word   | no     | binomial2
      n         | Binomial blending exponent      | scalar | no     | 2.0
    \endtable

    The inherited entries are elaborated in:
      - \link fixedValueFvPatchField.H \endlink
      - \link nutWallFunctionFvPatchScalarField.H \endlink

    Options for the \c blending entry:
    \verbatim
      stepwise    | Stepwise switch (discontinuous)
      max         | Maximum value switch (discontinuous)
      binomial2   | Binomial blending (smooth) n = 2
      binomial    | Binomial blending (smooth)
      exponential | Exponential blending (smooth)
      tanh        | Tanh blending (smooth)
    \endverbatim

    wherein \c omega predictions for the viscous and inertial sublayers are
    blended according to the following expressions:

    - \c stepwise:

    \f[
        \omega = \omega_{log} \qquad if \quad y^+ > y^+_{lam}
    \f]
    \f[
        \omega = \omega_{vis} \qquad if \quad y^+ <= y^+_{lam}
    \f]

    where
    \vartable
      \omega       | \f$\omega\f$ at \f$y^+\f$
      \omega_{vis} | \f$\omega\f$ computed by using viscous sublayer assumptions
      \omega_{log} |\f$\omega\f$ computed by using inertial sublayer assumptions
      y^+    | estimated wall-normal height of the cell centre in wall units
      y^+_{lam}  | estimated intersection of the viscous and inertial sublayers
    \endvartable


    - \c max (PH:Eq. 27):

    \f[
        \omega = max(\omega_{vis}, \omega_{log})
    \f]


    - \c binomial2 (ME:Eq. 15) (default):

    \f[
        \omega = \sqrt{(\omega_{vis})^2 + (\omega_{log})^2}
    \f]


    - \c binomial:

    \f[
        \omega = ((\omega_{vis})^n + (\omega_{log})^n)^{1/n}
    \f]

    where
    \vartable
        n             | Binomial blending exponent
    \endvartable


    - \c exponential (PH:Eq. 32):

    \f[
        \omega = \omega_{vis} \exp[-\Gamma] + \omega_{log} \exp[-1/\Gamma]
    \f]

    where (PH:Eq. 31)
    \vartable
        \Gamma | Blending expression
        \Gamma | \f$0.01 (y^+)^4 / (1.0 + 5.0 y^+)\f$
    \endvartable


    - \c tanh (KAS:Eqs. 33-34):

    \f[
        \omega = \phi \omega_{b1} + (1 - \phi)\omega_{b2}
    \f]

    where
    \vartable
        \phi        | \f$tanh((y^+/10)^4)\f$
        \omega_{b1} | \f$\omega_{vis} + \omega_{log}\f$
        \omega_{b2} | \f$(\omega_{vis}^{1.2} + \omega_{log}^1.2)^{1/1.2}\f$
    \endvartable


    \c G predictions for the viscous and inertial sublayers are blended
    in a stepwise manner, and \c G below \f$y^+_{lam}\f$ (i.e. in the viscous
    sublayer) is presumed to be zero.

Note
  - The coefficients \c Cmu, \c kappa, and \c E are obtained from
    the specified \c nutWallFunction in order to ensure that each patch
    possesses the same set of values for these coefficients.
  - The reason why \c binomial2 and \c binomial blending methods exist at
    the same time is to ensure the bitwise regression with the previous
    versions since \c binomial2 and \c binomial with \c n=2 will yield
    slightly different output due to the miniscule differences in the
    implementation of the basic functions (i.e. \c pow, \c sqrt, \c sqr).

See also
    - Foam::epsilonWallFunctionFvPatchScalarField

SourceFiles
    omegaWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef omegaWallFunctionFvPatchScalarField_H
#define omegaWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class turbulenceModel;

/*---------------------------------------------------------------------------*\
             Class omegaWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class omegaWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchField<scalar>
{
    // Private Enumerations

        //- Options for the blending treatment of viscous and inertial sublayers
        enum blendingType
        {
            STEPWISE,       //!< "Stepwise switch (discontinuous)"
            MAX,            //!< "Maximum value switch (discontinuous)"
            BINOMIAL2,      //!< "Binomial blending (smooth) n = 2"
            BINOMIAL,       //!< "Binomial blending (smooth)"
            EXPONENTIAL,    //!< "Exponential blending (smooth)"
            TANH            //!< "Tanh blending (smooth)"
        };

        //- Names for blendingType
        static const Enum<blendingType> blendingTypeNames;


    // Private Data

        //- Blending treatment (default = blendingType::BINOMIAL2)
        enum blendingType blending_;

        //- Blending exponent being used when
        //- blendingType is blendingType::BINOMIAL (default = 2)
        const scalar n_;


protected:

    // Protected Data

        //- Tolerance used in weighted calculations
        static scalar tolerance_;

        //- Deprecated(2019-11) Blending switch
        //  \deprecated(2019-11) - use blending:: options
        bool blended_;

        //- Initialised flag
        bool initialised_;

        //- Master patch ID
        label master_;

        //- beta1 coefficient
        scalar beta1_;

        //- Local copy of turbulence G field
        scalarField G_;

        //- Local copy of turbulence omega field
        scalarField omega_;

        //- List of averaging corner weights
        List<List<scalar>> cornerWeights_;


    // Protected Member Functions

        //- Set the master patch - master is responsible for updating all
        //- wall function patches
        virtual void setMaster();

        //- Create the averaging weights for cells which are bounded by
        //- multiple wall function faces
        virtual void createAveragingWeights();

        //- Helper function to return non-const access to an omega patch
        virtual omegaWallFunctionFvPatchScalarField& omegaPatch
        (
            const label patchi
        );

        //- Main driver to calculate the turbulence fields
        virtual void calculateTurbulenceFields
        (
            const turbulenceModel& turbulence,
            scalarField& G0,
            scalarField& omega0
        );

        //- Calculate the omega and G
        virtual void calculate
        (
            const turbulenceModel& turbulence,
            const List<scalar>& cornerWeights,
            const fvPatch& patch,
            scalarField& G,
            scalarField& omega
        );

        //- Return non-const access to the master patch ID
        virtual label& master()
        {
            return master_;
        }


public:

    //- Runtime type information
    TypeName("omegaWallFunction");


    // Constructors

        //- Construct from patch and internal field
        omegaWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        omegaWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- omegaWallFunctionFvPatchScalarField
        //- onto a new patch
        omegaWallFunctionFvPatchScalarField
        (
            const omegaWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        omegaWallFunctionFvPatchScalarField
        (
            const omegaWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new omegaWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        omegaWallFunctionFvPatchScalarField
        (
            const omegaWallFunctionFvPatchScalarField&,
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
                new omegaWallFunctionFvPatchScalarField(*this, iF)
            );
        }

    //- Destructor
    virtual ~omegaWallFunctionFvPatchScalarField() = default;


    // Member Functions

        // Access

            //- Return non-const access to the master's G field
            scalarField& G(bool init = false);

            //- Return non-const access to the master's omega field
            scalarField& omega(bool init = false);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Update the coefficients associated with the patch field
            virtual void updateWeightedCoeffs(const scalarField& weights);

            //- Manipulate matrix
            virtual void manipulateMatrix(fvMatrix<scalar>& matrix);

            //- Manipulate matrix with given weights
            virtual void manipulateMatrix
            (
                fvMatrix<scalar>& matrix,
                const scalarField& weights
            );


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
