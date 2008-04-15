/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2006-07 OpenCFD Ltd.
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
    Foam::massFlowRateInletVelocityFvPatchVectorField

Description
    Describes an massflow normal vector boundary condition by its magnitude
    as an integral over its area.
    The current density is used to correct the velocity.

    Example of the boundary condition specification:
    @verbatim
    inlet
    {
        type            massFlowRateInletVelocity;
        massFlowRate    0.2;             // Mass flow rate [kg/s]
        value           uniform (0 0 0); // placeholder
    }
    @endverbatim

Note
    - The value is positive inwards
    - May not work correctly for transonic inlets
    - Strange behaviour with potentialFoam since the U equation is not solved

SourceFiles
    massFlowRateInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef massFlowRateInletVelocityFvPatchVectorField_H
#define massFlowRateInletVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
               Class massFlowRateInletVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class massFlowRateInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Inlet integral mass flow rate
        scalar massFlowRate_;

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of the density field used to normalize the mass flux
        word rhoName_;


public:

   //- Runtime type information
   TypeName("massFlowRateInletVelocity");


   // Constructors

        //- Construct from patch and internal field
        massFlowRateInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        massFlowRateInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  massFlowRateInletVelocityFvPatchVectorField
        //  onto a new patch
        massFlowRateInletVelocityFvPatchVectorField
        (
            const massFlowRateInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        massFlowRateInletVelocityFvPatchVectorField
        (
            const massFlowRateInletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new massFlowRateInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        massFlowRateInletVelocityFvPatchVectorField
        (
            const massFlowRateInletVelocityFvPatchVectorField&,
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
                new massFlowRateInletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the mass flux
            scalar massFlowRate() const
            {
                return massFlowRate_;
            }

            //- Return reference to the mass flux to allow adjustment
            scalar& massFlowRate()
            {
                return massFlowRate_;
            }


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
