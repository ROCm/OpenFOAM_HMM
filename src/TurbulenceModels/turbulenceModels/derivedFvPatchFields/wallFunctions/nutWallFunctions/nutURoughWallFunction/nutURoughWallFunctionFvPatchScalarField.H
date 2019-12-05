/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::nutURoughWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a wall constraint on the turbulent
    kinematic viscosity, i.e. \c nut, when using wall functions for rough walls,
    based on velocity, \c U.

Usage
    \table
        Property     | Description             | Required    | Default value
        roughnessHeight | Roughness height     | yes         |
        roughnessConstant | Roughness constanr | yes         |
        roughnessFactor | Scaling factor       | yes         |
        Cmu        | Model coefficient         | no          | 0.09
        kappa      | Von Karman constant       | no          | 0.41
        E          | Model coefficient         | no          | 9.8
        maxIter    | Number of Newton-Raphson iterations  | no          | 10
        tolerance  | Convergence tolerance     | no          | 0.0001
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            nutURoughWallFunction;
        roughnessHeight     1e-5;
        roughnessConstant   0.5;
        roughnessFactor     1;

        // Optional entries
    }
    \endverbatim

Note
    Suffers from non-exact restart since correctNut() (called through
    turbulence->validate) returns a slightly different value every time
    it is called. See nutUSpaldingWallFunctionFvPatchScalarField.C.
    Can be avoided by seeding the NR with e.g. the laminar viscosity
    or tightening the convergence tolerance
    to e.g. 1e-7 and the max number of iterations to 100.

See also
    Foam::nutWallFunctionFvPatchScalarField

SourceFiles
    nutURoughWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutURoughWallFunctionFvPatchScalarField_H
#define nutURoughWallFunctionFvPatchScalarField_H

#include "nutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class nutURoughWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutURoughWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
    // Private Data

        // Roughness model parameters

            //- Height
            scalar roughnessHeight_;

            //- Constant
            scalar roughnessConstant_;

            //- Scale factor
            scalar roughnessFactor_;


        //- Max iterations in calcNut
        const label maxIter_;

        //- Convergence tolerance
        const scalar tolerance_;


    // Protected Member Functions

        //- Calculate yPlus
        virtual tmp<scalarField> calcYPlus(const scalarField& magUp) const;

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const;

        //- Write local wall function variables
        virtual void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("nutURoughWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutURoughWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutURoughWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- nutURoughWallFunctionFvPatchScalarField
        //- onto a new patch
        nutURoughWallFunctionFvPatchScalarField
        (
            const nutURoughWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutURoughWallFunctionFvPatchScalarField
        (
            const nutURoughWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nutURoughWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nutURoughWallFunctionFvPatchScalarField
        (
            const nutURoughWallFunctionFvPatchScalarField&,
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
                new nutURoughWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Return the roughness height
            scalar roughnessHeight() const
            {
                return roughnessHeight_;
            }

            //- Return reference to the roughness height to allow adjustment
            scalar& roughnessHeight()
            {
                return roughnessHeight_;
            }


            //- Return the roughness constant scale
            scalar roughnessConstant() const
            {
                return roughnessConstant_;
            }

            //- Return reference to the roughness constant to allow adjustment
            scalar& roughnessConstant()
            {
                return roughnessConstant_;
            }

            //- Return the roughness scale factor
            scalar roughnessFactor() const
            {
                return roughnessFactor_;
            }

            //- Return reference to the roughness scale factor to allow
            //- adjustment
            scalar& roughnessFactor()
            {
                return roughnessFactor_;
            }


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
