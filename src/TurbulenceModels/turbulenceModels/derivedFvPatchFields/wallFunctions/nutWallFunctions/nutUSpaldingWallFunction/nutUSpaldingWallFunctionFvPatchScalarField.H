/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::nutUSpaldingWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a turbulent kinematic viscosity condition
    when using wall functions for rough walls, based on velocity,  using
    Spalding's law to give a continuous nut profile to the wall (y+ = 0)

        \f[
            y^+ = u^+ + \frac{1}{E} \left[exp(\kappa u^+) - 1 - \kappa u^+\,
                - 0.5 (\kappa u^+)^2 - \frac{1}{6} (\kappa u^+)^3\right]
        \f]

    where
    \vartable
        y^+     | non-dimensional position
        u^+     | non-dimensional velocity
        \kappa  | Von Karman constant
    \endvartable


Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            nutUSpaldingWallFunction;

        // .. nutWallFunctionFvPatchScalarField entries ..

    }
    \endverbatim

    Additional entries (see below)
    \table
        Property  | Description               | Required   | Default value
        maxIter   | Number of N-R iterations  | no         | 10
        tolerance | Convergence tolerance     | no         | 0.01
    \endtable

See also
    Foam::nutWallFunctionFvPatchScalarField

Note
    Suffers from non-exact restart since correctNut() (called through
    turbulence->validate) returns a slightly different value every time
    it is called. This is since the seed for the Newton-Raphson iteration
    uses the current value of *this (= nut). Can be avoided by seeding the
    NR with e.g. the laminar viscosity or tightening the convergence tolerance
    to e.g. 1e-7 and the max number of iterations to 100.

SourceFiles
    nutUSpaldingWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutUSpaldingWallFunctionFvPatchScalarField_H
#define nutUSpaldingWallFunctionFvPatchScalarField_H

#include "nutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class nutUSpaldingWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutUSpaldingWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
protected:

    // Protected data

        //- Max iterations in calcNut
        const label maxIter_;

        //- Convergence tolerance
        const scalar tolerance_;

        //- Uncomment in case of intrumentation
        //mutable uint64_t invocations_;
        //mutable uint64_t nontrivial_;
        //mutable uint64_t nonconvergence_;
        //mutable uint64_t iterations_;


    // Protected Member Functions

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const;

        //- Calculate the friction velocity
        virtual tmp<scalarField> calcUTau(const scalarField& magGradU) const;

        //- Write local wall function variables
        virtual void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("nutUSpaldingWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutUSpaldingWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutUSpaldingWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutUSpaldingWallFunctionFvPatchScalarField
        //  onto a new patch
        nutUSpaldingWallFunctionFvPatchScalarField
        (
            const nutUSpaldingWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutUSpaldingWallFunctionFvPatchScalarField
        (
            const nutUSpaldingWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nutUSpaldingWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nutUSpaldingWallFunctionFvPatchScalarField
        (
            const nutUSpaldingWallFunctionFvPatchScalarField&,
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
                new nutUSpaldingWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    //- Destructor
    virtual ~nutUSpaldingWallFunctionFvPatchScalarField();


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
