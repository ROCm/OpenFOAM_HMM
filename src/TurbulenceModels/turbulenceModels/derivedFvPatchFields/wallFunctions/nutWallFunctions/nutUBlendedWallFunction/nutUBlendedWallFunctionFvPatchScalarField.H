/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2019 OpenCFD Ltd.
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
    Foam::nutUBlendedWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a turbulent kinematic viscosity condition
    when using wall functions, based on a blending of laminar sub-layer and
    log region contributions.

        \f[
            u_\tau = (u_{\tau,v}^n + u_{\tau,l}^n)^{1/n}
        \f]

    where
    \vartable
        u_\tau     | friction velocity
        u_{\tau,v} | friction velocity in the viscous region
        u_{\tau,l} | friction velocity in the log region
    \endvartable


Usage

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            nutUBlendedWallFunction;
    }
    \endverbatim

    Reference:
    See the section that describes 'automatic wall treatment'
    \verbatim
        Menter, F., Carregal Ferreira, J., Esch, T., Konno, B. (2003).
        The SST Turbulence Model with Improved Wall Treatment
        for Heat Transfer Predictions in Gas Turbines.
        Proceedings of the International Gas Turbine Congress 2003 Tokyo
    \endverbatim

Note
    The full 'automatic wall treatment' description also requires use of the
    Foam::omegaWallFunction with the \c blended flag set to 'on'

    Suffers from non-exact restart since correctNut() (called through
    turbulence->validate) returns a slightly different value every time
    it is called. See nutUSpaldingWallFunctionFvPatchScalarField.C

SeeAlso
    Foam::nutWallFunctionFvPatchScalarField
    Foam::omegaWallFunctionFvPatchScalarField

SourceFiles
    nutUBlendedWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutUBlendedWallFunctionFvPatchScalarField_H
#define nutUBlendedWallFunctionFvPatchScalarField_H

#include "nutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class nutUBlendedWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutUBlendedWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
protected:

    // Protected data

        //- Model coefficient; default = 4
        scalar n_;


    // Protected Member Functions

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const;

        //- Calculate the friction velocity
        virtual tmp<scalarField> calcUTau(const scalarField& magGradU) const;


public:

    //- Runtime type information
    TypeName("nutUBlendedWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutUBlendedWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutUBlendedWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutUBlendedWallFunctionFvPatchScalarField
        //  onto a new patch
        nutUBlendedWallFunctionFvPatchScalarField
        (
            const nutUBlendedWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutUBlendedWallFunctionFvPatchScalarField
        (
            const nutUBlendedWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nutUBlendedWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nutUBlendedWallFunctionFvPatchScalarField
        (
            const nutUBlendedWallFunctionFvPatchScalarField&,
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
                new nutUBlendedWallFunctionFvPatchScalarField(*this, iF)
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

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
