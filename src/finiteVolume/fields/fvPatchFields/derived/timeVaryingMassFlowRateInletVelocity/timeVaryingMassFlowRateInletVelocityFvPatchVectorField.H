/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2006-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::timeVaryingMassFlowRateInletVelocityFvPatchVectorField

Description
    A time-varying form of a massflow normal vector boundary condition.

    Example of the boundary condition specification:
    @verbatim
    inlet
    {
        type            timeVaryingMassFlowRateInletVelocity;
        massFlowRate    0.2;              // Massflow rate [kg/s]
        value           uniform (0 0 0);  // placeholder
        timeDataFile   "time-series";
        timeBounding    repeat;           // (error|warn|clamp|repeat)
    }
    @endverbatim

Note
    - The value is positive inwards
    - may not work correctly for transonic inlets!
    - strange behaviour with potentialFoam since the U equation is not solved

See Also
    Foam::timeSeries and Foam::massFlowRateInletVelocityFvPatchVectorField

SourceFiles
    timeVaryingMassFlowRateInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingMassFlowRateInletVelocityFvPatchVectorField_H
#define timeVaryingMassFlowRateInletVelocityFvPatchVectorField_H

#include "massFlowRateInletVelocityFvPatchVectorField.H"
#include "timeSeries.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
        Class timeVaryingMassFlowRateInletVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingMassFlowRateInletVelocityFvPatchVectorField
:
    public massFlowRateInletVelocityFvPatchVectorField
{
    // Private data

        //- file containing time/massFlowRate
        fileName timeDataFile_;

        //- the time series being used, including the bounding treatment
        timeSeries<scalar> timeSeries_;

        //- interpolate the value at the current time
        scalar currentValue();

public:

    //- Runtime type information
    TypeName("timeVaryingMassFlowRateInletVelocity");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingMassFlowRateInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingMassFlowRateInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch field onto a new patch
        timeVaryingMassFlowRateInletVelocityFvPatchVectorField
        (
            const timeVaryingMassFlowRateInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        timeVaryingMassFlowRateInletVelocityFvPatchVectorField
        (
            const timeVaryingMassFlowRateInletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new timeVaryingMassFlowRateInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        timeVaryingMassFlowRateInletVelocityFvPatchVectorField
        (
            const timeVaryingMassFlowRateInletVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new timeVaryingMassFlowRateInletVelocityFvPatchVectorField(*this, iF)
            );
        }

    // Member functions

        // Access

            //- Return the out-of-bounds treatment as a word
            word timeBounding() const
            {
                return timeSeries_.bounding();
            }

            //- Return the time series used
            const timeSeries<scalar>& timeData() const
            {
                return timeSeries_;
            }

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
