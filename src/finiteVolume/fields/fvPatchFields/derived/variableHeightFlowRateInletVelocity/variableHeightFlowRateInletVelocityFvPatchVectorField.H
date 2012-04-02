/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::variableHeightFlowRateInletVelocityFvPatchVectorField

Description
    Describes a volumetric/mass flow normal vector boundary condition by its
    magnitude as an integral over its area.

    The basis of the patch (volumetric or mass) is determined by the
    dimensions of the flux, phi.
    The current density is used to correct the velocity when applying the
    mass basis.

    The flow rate is made proportional to the phase fraction alpha at each face
    of the patch and alpha is ensured to be bound between 0 and 1.

    Example of the boundary condition specification:
    \verbatim
    inlet
    {
        type            variableHeightFlowRateInletVelocity;
        flowRate        0.2;        // Volumetric/mass flow rate [m3/s or kg/s]
        value           uniform (0 0 0); // placeholder
    }
    \endverbatim

Note
    - The value is positive inwards
    - May not work correctly for transonic inlets
    - Strange behaviour with potentialFoam since the U equation is not solved

SourceFiles
    variableHeightFlowRateInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef variableHeightFlowRateInletVelocityFvPatchVectorField_H
#define variableHeightFlowRateInletVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
          Class variableHeightFlowRateInletVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class variableHeightFlowRateInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Inlet integral flow rate
        scalar flowRate_;

public:

   //- Runtime type information
   TypeName("variableHeightFlowRateInletVelocity");


   // Constructors

        //- Construct from patch and internal field
        variableHeightFlowRateInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        variableHeightFlowRateInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  variableHeightFlowRateInletVelocityFvPatchVectorField
        //  onto a new patch
        variableHeightFlowRateInletVelocityFvPatchVectorField
        (
            const variableHeightFlowRateInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        variableHeightFlowRateInletVelocityFvPatchVectorField
        (
            const variableHeightFlowRateInletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new variableHeightFlowRateInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        variableHeightFlowRateInletVelocityFvPatchVectorField
        (
            const variableHeightFlowRateInletVelocityFvPatchVectorField&,
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
                new variableHeightFlowRateInletVelocityFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Access

            //- Return the flux
            scalar flowRate() const
            {
                return flowRate_;
            }

            //- Return reference to the flux to allow adjustment
            scalar& flowRate()
            {
                return flowRate_;
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
