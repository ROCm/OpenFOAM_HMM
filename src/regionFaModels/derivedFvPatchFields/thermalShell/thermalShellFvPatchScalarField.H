/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::compressible::thermalShellFvPatchScalarField

Group
    grpThermoBoundaryConditions

Description
    This boundary condition provides a coupled temperature condition between
    a primary region (3D mesh) and a thermal shell model (2D mesh).

    The primary region boundary creates the finite area region
    and evolves its energy equation.

Usage
    Example of the boundary condition specification:
    \verbatim
    <masterPatchName>
    {
        // Mandatory entries (unmodifiable)
        type                compressible::thermalShell;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property | Description                           | Type  | Reqd | Dflt
      type     | Type name: compressible::thermalShell | word  | yes  | -
    \endtable

    The inherited entries are elaborated in:
      - \link fixedValueFvPatchField.H \endlink
      - \link thermalShellModel.H \endlink

Note
  - The two-dimensional area mesh needs to be
    generated in the pre-processing steps.

SourceFiles
    thermalShellFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef thermalShellFvPatchScalarField_H
#define thermalShellFvPatchScalarField_H

#include "autoPtr.H"
#include "thermalShellModel.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
             Class thermalShellFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class thermalShellFvPatchScalarField
:
    public fixedValueFvPatchField<scalar>
{
    // Private Data

        //- Thermal baffle
        autoPtr<regionModels::thermalShellModel> baffle_;

        //- Dictionary
        dictionary dict_;


public:

    //- Runtime type information
    TypeName("compressible::thermalShell");


    // Constructors

        //- Construct from patch and internal field
        thermalShellFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        thermalShellFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- thermalShellFvPatchScalarField onto a new patch
        thermalShellFvPatchScalarField
        (
            const thermalShellFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new thermalShellFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        thermalShellFvPatchScalarField
        (
            const thermalShellFvPatchScalarField&,
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
                new thermalShellFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
