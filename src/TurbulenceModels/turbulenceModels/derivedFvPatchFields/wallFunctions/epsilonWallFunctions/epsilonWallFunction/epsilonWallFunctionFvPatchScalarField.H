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
    Foam::epsilonWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a wall constraint on the turbulent kinetic
    energy dissipation rate, i.e. \c epsilon, and the turbulent kinetic
    energy production contribution, i.e. \c G, for low- and high-Reynolds
    number turbulence models.

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
    \endverbatim

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries (unmodifiable)
        type             epsilonWallFunction;

        // Optional entries (unmodifiable)
        lowReCorrection  false;
        blending         stepwise;
        n                2.0;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                       | Type | Req'd  | Dflt
      type      | Type name: epsilonWallFunction    | word | yes    | -
      lowReCorrection | Flag: apply low-Re correction | bool |  no  | false
      blending  | Viscous/inertial sublayer blending method <!--
                                                --> | word |  no    | stepwise
      n         | Binomial blending exponent        | scalar | no   | 2.0
    \endtable

    The inherited entries are elaborated in:
      - \link fixedValueFvPatchField.H \endlink
      - \link nutWallFunctionFvPatchScalarField.H \endlink

    Options for the \c blending entry:
    \verbatim
      stepwise    | Stepwise switch (discontinuous)
      max         | Maximum value switch (discontinuous)
      binomial    | Binomial blending (smooth)
      exponential | Exponential blending (smooth)
    \endverbatim

    wherein \c epsilon predictions for the viscous and inertial sublayers are
    blended according to the following expressions:

    - \c stepwise (default):

    \f[
        \epsilon = \epsilon_{vis} \qquad if \quad y^+ < y^+_{lam}
    \f]
    \f[
        \epsilon = \epsilon_{log} \qquad if \quad y^+ >= y^+_{lam}
    \f]

    where
    \vartable
      \epsilon       | \f$\epsilon\f$ at \f$y^+\f$
      \epsilon_{vis} | \f$\epsilon\f$ computed by viscous subl. assumptions
      \epsilon_{log} | \f$\epsilon\f$ computed by inertial subl. assumptions
      y^+    | estimated wall-normal height of the cell centre in wall units
      y^+_{lam} | estimated intersection of the viscous and inertial sublayers
    \endvartable


    - \c max (PH:Eq. 27):

    \f[
        \epsilon = max(\epsilon_{vis}, \epsilon_{log})
    \f]


    - \c binomial (ME:Eqs. 15-16):

    \f[
        \epsilon = ((\epsilon_{vis})^n + (\epsilon_{log})^n)^{1/n}
    \f]
    where
    \vartable
      n               | Binomial blending exponent
    \endvartable


    - \c exponential (PH:Eq. 32):

    \f[
        \epsilon = \epsilon_{vis} \exp[-\Gamma] +\epsilon_{log} \exp[-1/\Gamma]
    \f]
    where (PH:p. 193)
    \vartable
      \Gamma_\epsilon | \f$\Gamma = 0.001 (y^+)^4 / (1.0 + y^+)\f$
      \Gamma_G        | \f$\Gamma = 0.01  (y^+)^4 / (1.0 + 5.0 y^+)\f$
      \Gamma_\epsilon | Blending expression for \f$\epsilon\f$
      \Gamma_G        | Blending expression for \f$G\f$
    \endvartable


    \c G predictions for the viscous and inertial sublayers are blended
    in a stepwise manner, and \c G below \f$y^+_{lam}\f$ (i.e. in the viscous
    sublayer) is presumed to be zero.

Note
  - The coefficients \c Cmu, \c kappa, and \c E are obtained from
    the specified \c nutWallFunction in order to ensure that each patch
    possesses the same set of values for these coefficients.
  - \c lowReCorrection operates with only \c stepwise blending treatment to
    ensure the backward compatibility.
  - If \c lowReCorrection is \c on, \c stepwise blending treatment is fully
    active.
  - If \c lowReCorrection is \c off, only the inertial sublayer prediction
    is used in the wall function, hence high-Re mode operation.

See also
    - Foam::omegaWallFunctionFvPatchScalarField

SourceFiles
    epsilonWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef epsilonWallFunctionFvPatchScalarField_H
#define epsilonWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class turbulenceModel;

/*---------------------------------------------------------------------------*\
           Class epsilonWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class epsilonWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchField<scalar>
{
    // Private Enumerations

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


    // Private Data

        //- Blending treatment (default = blendingType::STEPWISE)
        const enum blendingType blending_;

        //- Binomial blending exponent being used when
        //- blendingType is blendingType::BINOMIAL (default = 2)
        const scalar n_;


protected:

    // Protected Data

        //- Tolerance used in weighted calculations
        static scalar tolerance_;

        //- Apply low-Re correction term (default = no)
        const bool lowReCorrection_;

        //- Initialised flag
        bool initialised_;

        //- Master patch ID
        label master_;

        //- Local copy of turbulence G field
        scalarField G_;

        //- Local copy of turbulence epsilon field
        scalarField epsilon_;

        //- List of averaging corner weights
        List<List<scalar>> cornerWeights_;


    // Protected Member Functions

        //- Set the master patch - master is responsible for updating all
        //- wall function patches
        virtual void setMaster();

        //- Create the averaging weights for cells which are bounded by
        //- multiple wall function faces
        virtual void createAveragingWeights();

        //- Helper function to return non-const access to an epsilon patch
        virtual epsilonWallFunctionFvPatchScalarField& epsilonPatch
        (
            const label patchi
        );

        //- Main driver to calculate the turbulence fields
        virtual void calculateTurbulenceFields
        (
            const turbulenceModel& turbulence,
            scalarField& G0,
            scalarField& epsilon0
        );

        //- Calculate the epsilon and G
        virtual void calculate
        (
            const turbulenceModel& turbulence,
            const List<scalar>& cornerWeights,
            const fvPatch& patch,
            scalarField& G,
            scalarField& epsilon
        );

        //- Return non-const access to the master patch ID
        virtual label& master()
        {
            return master_;
        }


public:

    //- Runtime type information
    TypeName("epsilonWallFunction");


    // Constructors

        //- Construct from patch and internal field
        epsilonWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        epsilonWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- epsilonWallFunctionFvPatchScalarField
        //- onto a new patch
        epsilonWallFunctionFvPatchScalarField
        (
            const epsilonWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        epsilonWallFunctionFvPatchScalarField
        (
            const epsilonWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new epsilonWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        epsilonWallFunctionFvPatchScalarField
        (
            const epsilonWallFunctionFvPatchScalarField&,
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
                new epsilonWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    //- Destructor
    virtual ~epsilonWallFunctionFvPatchScalarField() = default;


    // Member Functions

        // Access

            //- Return non-const access to the master's G field
            scalarField& G(bool init = false);

            //- Return non-const access to the master's epsilon field
            scalarField& epsilon(bool init = false);


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
